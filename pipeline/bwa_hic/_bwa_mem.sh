
ref=$1
ref_base=`basename $ref`
ref_base=${ref_base/.fasta/}
ref_base=${ref_base/.fa/}
ref_base=${ref_base/.fna/}

threads=$2

export BWA=/data/projects/phillippy/software/bwa/bwa-0.7.15
export PATH=$BWA:$PATH

if [ ! -e $ref.sa ] ; then
        echo "No bwa index file found. Indexing..."
        bwa index $ref
        echo "Done"
        echo ""
fi

for i in 3 5 2     #seq 2 5
do
	./_sge.sh $threads bwa_mem_$i ./_bwa_mem_parallel.sh "$ref $i $threads"
done
