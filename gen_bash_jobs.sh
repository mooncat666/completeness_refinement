PIPE_DIR=/g/scb2/bork/groudko/scripts/ref_v1.1
SCRIPT=${PIPE_DIR}/jobs.sh

rm $SCRIPT
touch $SCRIPT

for s in "'Fusobacterium nucleatum'" "'Fusobacterium nucleatum_J'" "'ANI_95_v1_000547972'" "'Clostridioides difficile'" "'ANI_95_v1_000537534'"
do
# "'ANI_95_v1_000537534'"
echo "bash ${PIPE_DIR}/pipe.sh -s $s -d /g/scb/bork/groudko/ref_v1.1 2>>pipe.log &" >> $SCRIPT
done

echo "disown -a" >> $SCRIPT


