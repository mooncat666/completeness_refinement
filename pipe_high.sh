# usage: bash pipe.sh -s species_name [-d outdir --temp temp_dir -h]
# add a scratch location for all intermediate files

path=/g/scb2/bork/groudko/scripts/ref_v1.1 # path to scripts
source ~/.bashrc
outdir='.'
final='.'

# set options

while [ -n "$1" ] ; do
case "$1" in
-d) outdir=$2 
final=$2 ;;
-s | --species) species=$2 ;;
--temp) temp_outdir=$2 ;;
-h) echo "Usage: bash -i pipe.sh -s species [-d output_directory] [--temp temp_dir] [-h]"
exit 1
shift ;;
esac
shift
done

if [[ -n ${temp_outdir} ]]; # --temp option is set
then
outdir=${temp_outdir}
fi

if [[ -z ${species} ]];
then
echo "Mandatory option -s (species name or ANI cluster ID) is not set"
exit 1
fi

if [[ ! -d $outdir ]];
then
mkdir $outdir
fi


# 0 get data
 
cd $outdir

mkdir 0_data
module purge
conda activate jupyternotebook 
# srun --mem 32000 -c 1 -t 600 --qos highest python ${path}/0_get_data.py "${species}"
conda deactivate

species=$( echo $species | sed 's/ /_/g' )


mkdir 0_data/${species}

if [[ ! -d 0_data/${species}/unzipped_init_MAGs ]];
then

mkdir 0_data/${species}/unzipped_init_MAGs
mkdir 0_data/${species}/unzipped_samples

# gunzip files

for f in 0_data/${species}/init_MAGs/*fa.gz
do
gunzip -c $f > 0_data/${species}/unzipped_init_MAGs/$(basename $f .gz)
done

for f in 0_data/${species}/samples/*fa.gz
do
gunzip -c $f > 0_data/${species}/unzipped_samples/$(basename $f .gz)
done

fi


mkdir 0_data/${species}/gunc
mkdir 0_data/${species}/checkm2

# 1.1 run gunc refinement

module load gunc

mkdir 1.1_gunc_refined
mkdir 1.1_gunc_refined/${species}
mkdir 1.1_gunc_refined/${species}/gunc

# run gunc
srun -n 1 -c 16 -t 2000 --mem 64000 gunc run -d 0_data/${species}/unzipped_init_MAGs -o 1.1_gunc_refined/${species}/gunc -e .fa -t 16 
cp  1.1_gunc_refined/${species}/gunc/GUNC.maxCSS_level.tsv 0_data/${species}/gunc/GUNC.maxCSS_level.tsv # for 4.2
cd 1.1_gunc_refined/${species}/gunc/diamond_output
for f in *diamond.out ; do f=$(basename $f .diamond.out) ; awk -v f=$f '{$1=f":"$1 ; print $0}' $f.diamond.out > temp && mv temp $f.diamond.out ; done
cd ../gene_calls
for file in *.faa ; do sed -e "s/>/>$(basename -s .faa ${file}):/g" $file > temp && mv temp $file ; done
cd ../../../..

# run refinement
srun -n 1 -c 8 --mem 32000 -t 1000 python /g/scb2/bork/groudko/scripts/refine_gunc.py 1.1_gunc_refined/${species}/gunc
python /g/scb2/bork/groudko/scripts/step2_refine_gunc.py 0_data/${species}/unzipped_init_MAGs 1.1_gunc_refined/${species}/gunc
gunc_clade=$( basename $( ls -d 1.1_gunc_refined/${species}/gunc/refinement/refined_genomes/*_*/ ))

mkdir 1.1_gunc_refined/${species}/gunc/refinement/refined_genomes/${gunc_clade}/gunc
# for 4.2
srun -n 1 -c 16 -t 2000 --mem 64000 gunc run -d 1.1_gunc_refined/${species}/gunc/refinement/refined_genomes/${gunc_clade} -o 1.1_gunc_refined/${species}/gunc/refinement/refined_genomes/${gunc_clade}/gunc -e .fa -t 16 &

# run checkm to estimate completeness

module purge
module load checkm2

# on initial data for 4.2
srun -n 1 -c 32 --mem 64000 -t 1200 checkm2 predict -i 0_data/${species}/unzipped_init_MAGs -t 32 -o 0_data/${species}/checkm2 -x .fa &

srun -n 1 -c 32 --mem 64000 -t 1200 checkm2 predict -i 1.1_gunc_refined/${species}/gunc/refinement/refined_genomes/${gunc_clade} -t 32 -o 1.1_gunc_refined/${species}/gunc/refinement/refined_genomes/${gunc_clade}/checkm2 -x .fa 
module purge


gunc_clade=$( basename $( ls -d 1.1_gunc_refined/${species}/gunc/refinement/refined_genomes/*_*/ ))


conda activate jupyternotebook
python ${path}/add_gunc_results.py 1.1_gunc_refined/${species}/gunc/refinement/refined_genomes/${gunc_clade}
conda deactivate

mkdir 1.1_gunc_refined/${species}/compl50
for f in $( awk -F'\t' '{ if ($2>=50 && $3<5 && $4=="True" && NR>1) {print $1} }' 1.1_gunc_refined/${species}/gunc/refinement/refined_genomes/${gunc_clade}/checkm2/quality_report_gunc.tsv )
do
ln -s -r 1.1_gunc_refined/${species}/gunc/refinement/refined_genomes/${gunc_clade}/$f*.fa 1.1_gunc_refined/${species}/compl50
done

# check if successful
# implement: check return codes
if [ ! "$(ls -A 1.1_gunc_refined/${species}/compl50)" ]; then
echo "ERROR: no refined bins are passing filtering. Terminating script execution on species ${species}" >> /dev/stderr
exit 1
fi

# 1.2 deduplicate

module purge
module load BBMap

mkdir 1.2_dedupl
mkdir 1.2_dedupl/${species}

cat 1.1_gunc_refined/${species}/compl50/*fa > 1.2_dedupl/${species}/merged.fa
srun --mem 64000 -n 1 -c 16 -t 1200 dedupe.sh threads=16 in=1.2_dedupl/${species}/merged.fa out=1.2_dedupl/${species}/dd_merged.fa minidentity=99 exact=f


# 1.3 discard short contigs from assemblies

mkdir 1.3_discard_short
mkdir 1.3_discard_short/${species}

for f in 0_data/${species}/unzipped_samples/*fa
do
reformat.sh in=$f out=1.3_discard_short/${species}/$(basename $f) minlength=1500 2>/dev/null
done

module purge



# 2.1 build variation graph

mkdir 2.1_build_vg
mkdir 2.1_build_vg/${species}

conda activate pggb
module load SAMtools
samtools faidx 1.2_dedupl/${species}/dd_merged.fa
module purge
files_number=$( ls 1.1_gunc_refined/${species}/compl50/*fa | wc -l )

srun -n 1 -c 16 -t 2000 --mem 200000 pggb -i 1.2_dedupl/${species}/dd_merged.fa -o 2.1_build_vg/${species} -P asm20 -t 16 -p 94 -s 10000 -C cons,100 -n ${files_number} -k 49 -G 7919,8069 --skip-viz || check=True
# rerun if failed 
# check $check
# check if consensus ???
conda deactivate


#rm with comment
#gunc_clade=$( basename $( ls -d 1.1_gunc_refined/${species}/gunc/refinement/refined_genomes/*_*/ ))



# 3.1 & 3.2 map contigs from assemblies to graph

mkdir 3.1_map_contigs
mkdir 3.1_map_contigs/${species}
mkdir 3.2_filter_alignments
mkdir 3.2_filter_alignments/${species}

# 3.1 map contigs & 3.2 filter alignments

bash ${path}/gen_jobs.sh $species $PWD
module load submitjobs
submitjob -s SLURM -n ${species:0:4} -c 4 -m 64 -a 100 -k 48:00:00 3.1_map_contigs/${species}/jobs_map.sh | (read job_id ; python ${path}/check_job_status.py $job_id ; echo $job_id > 3.1_map_contigs/${species}/job_id.txt )

# rerun failed
rm 3.1_map_contigs/${species}/rerun_jobs_map.sh
touch 3.1_map_contigs/${species}/rerun_jobs_map.sh
sacct -j $( cat 3.1_map_contigs/${species}/job_id.txt ) -s nf,f,ca --format=jobid%20 | sed 1,2d | cut -d "_" -f 2 | sed -r "s/\..*//" | tr -d ' ' | sort -u | while read n ; do sed -n "$np" 3.1_map_contigs/${species}/jobs_map.sh >> 3.1_map_contigs/${species}/rerun_jobs_map.sh ; done
if [ -s 3.1_map_contigs/${species}/rerun_jobs_map.sh ]; then
submitjob -s SLURM -n ${species:0:4} -c 4 -m 64 -a 100 -k 48:00:00 3.1_map_contigs/${species}/rerun_jobs_map.sh | (read job_id ; python ${path}/check_job_status.py $job_id)
fi
# rerun timeout
# sacct -j $( cat 3.1_map_contigs/${species}/job_id.txt ) -s to --format=jobid%20 | sed 1,2d | cut -d "_" -f 2 | sed -r "s/\..*//" | tr -d ' ' | sort -u

module purge



# 4.1 gunc refine mapped bins

module load gunc

mkdir 4.1_refine_output
mkdir 4.1_refine_output/${species}
mkdir 4.1_refine_output/${species}/gunc

# run gunc
srun -n 1 -c 16 -t 2000 --mem 64000 gunc run -d 3.2_filter_alignments/${species} -o 4.1_refine_output/${species}/gunc -e .fa -t 16
cd 4.1_refine_output/${species}/gunc/diamond_output
for f in *diamond.out ; do f=$(basename $f .diamond.out) ; awk -v f=$f '{$1=f":"$1 ; print $0}' $f.diamond.out > temp && mv temp $f.diamond.out ; done
cd ../gene_calls
for file in *.faa ; do sed -e "s/>/>$(basename -s .faa ${file}):/g" $file > temp && mv temp $file ; done
cd ../../../..
# run refinement
srun -n 1 -c 8 --mem 32000 -t 600 python /g/scb2/bork/groudko/scripts/refine_gunc.py 4.1_refine_output/${species}/gunc
python /g/scb2/bork/groudko/scripts/step2_refine_gunc.py 3.2_filter_alignments/${species} 4.1_refine_output/${species}/gunc
gunc_clade=$( basename $( ls -d 4.1_refine_output/${species}/gunc/refinement/refined_genomes/*_*/ ))

# 4.2 create stats

# run checkm and gunc on refined output and initial bins

mkdir 4.2_output_n_stats
mkdir 4.2_output_n_stats/${species}
mkdir 4.2_output_n_stats/${species}/gunc
mkdir 4.2_output_n_stats/${species}/checkm2

for f in 4.1_refine_output/${species}/gunc/refinement/refined_genomes/${gunc_clade}/*fa
do
ln -s -r $f 4.2_output_n_stats/${species}/$(basename $f)
done

cd 4.2_output_n_stats/${species}

module purge
module load gunc
srun -n 1 -c 16 -t 2000 --mem 64000 gunc run -d . -o gunc -e .fa -t 16
module purge

module load checkm2
srun -n 1 -c 32 --mem 64000 -t 1200 checkm2 predict -i .  -t 32 -o checkm2 -x .fa
module purge

# create plots
conda activate jupyternotebook
python ${path}/merge_gunc_checkm2_3.py ../../0_data/${species} 1.1_gunc_refined/${species}/gunc/refinement/refined_genomes/${gunc_clade}
Rscript ${path}/plots.r
conda deactivate
module purge

# write end files and stats to FINAL directory

if [ $final != $outdir ]; then
mkdir ${final}/4.2_output_n_stats
mkdir ${final}/4.2_output_n_stats/${species}
cp *pdf *tsv ${final}/4.2_output_n_stats/${species}
fi

echo "end"
