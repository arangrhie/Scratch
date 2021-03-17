DIPLOID="--diploid"
ALGORITHM=arrow
prefix=$1
asm=$prefix.fasta	# fasta
jobid=$SGE_TASK_ID
chunk_num=`echo $jobid | awk '{print ($1-1)}'`
chunk=$prefix/by_contigs/$prefix.chunk$chunk_num.xml
pe_slots=12

mkdir -p $prefix/out

echo "=== Start running variantCaller --algorithm=arrow on $chunk ==="
echo "variantCaller --skipUnrecognizedContigs $DIPLOID -x 5 -q 20 -X120 --log-level=ERROR -j$pe_slots --algorithm=$ALGORITHM -r $asm -o $prefix/out/$jobid.gff -o $prefix/out/$jobid.fastq -o $prefix/out/$jobid.fasta $chunk"
variantCaller --skipUnrecognizedContigs $DIPLOID -x 5 -q 20 -X120 --log-level=ERROR -j$pe_slots --algorithm=$ALGORITHM -r $asm -o $prefix/out/$jobid.gff -o $prefix/out/$jobid.fastq -o $prefix/out/$jobid.fasta $chunk
echo ""
