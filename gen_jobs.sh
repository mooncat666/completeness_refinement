# Usage: ./gen_jobs.sh species pipe_outdir
SPECIES=$1
OUTDIR=$2
CONS=$3
SCRIPT=${OUTDIR}/3.1_map_contigs/${SPECIES}/jobs_map.sh
rm $SCRIPT
touch $SCRIPT



for FILE in ${OUTDIR}/1.3_discard_short/${SPECIES}/*fa
do
echo "
bash /g/scb2/bork/groudko/scripts/ref_v1.1/map_contigs.sh -d $OUTDIR -s $SPECIES -c $CONS -f $FILE 
" >> $SCRIPT
done
