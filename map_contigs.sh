script_path=/g/scb2/bork/groudko/scripts/ref_v1.1


# set options

while [ -n "$1" ] ; do
case "$1" in
-d) outdir=$2 ;;
-s | --species) species=$2 ;;
-c | --consensus) cons=$2;;
-f) file=$2;;
-h) echo "Usage: bash map_contigs.sh -d outdir -s species -c consensus_block_length"
exit 1
shift ;;
esac
shift
done


graph=${outdir}/2.2_edit_vg/${species}/filtered_graph.gfa

init_bin_dir=$( ls -d ${outdir}/1.1_gunc_refined/${species}/gunc/refinement/refined_genomes/*_*/ )
init_bin=$( ls ${init_bin_dir} | grep -m 1 $(basename $file -assembled.fa) )

source ~/.bashrc


# map contigs from each sample to the graph

conda activate graphaligner_clone
GraphAligner --min-alignment-score 100 --multimap-score-fraction 0 -a ${outdir}/3.1_map_contigs/${species}/mapped_$(basename $file .fa).gaf -x vg -f $file -g ${graph} -t 4 --precise-clipping 0.95
conda deactivate


# filter mappings

conda activate jupyternotebook_clone
Rscript ${script_path}/filter_mappings.r ${outdir}/3.1_map_contigs/${species}/mapped_$(basename $file .fa).gaf 0.9 | python ${script_path}/add_init.py -i ${init_bin_dir}${init_bin} -s $file | python /g/scb2/bork/groudko/scripts/fasta_from_ids.py $file > ${outdir}/3.2_filter_alignments/${species}/mapped_$(basename $file )
conda deactivate
