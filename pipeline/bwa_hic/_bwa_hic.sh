ref=$1
ref_base=`basename $ref`
ref_base=${ref_base/.fasta/}
ref_base=${ref_base/.fa/}
ref_base=${ref_base/.fna/}

export BWA=/data/projects/phillippy/software/bwa/bwa-0.7.15
export PATH=$BWA:$PATH

# for checking version we are running
bwa
samtools


if [ ! -e $ref.sa ] ; then
	echo "No bwa index file found. Indexing..."
	bwa index $ref
	echo "Done"
	echo ""
fi

threads=30
in_map=hic_reads.map

for i in $(seq 2 5)	#seq 2 5
do
	tech=`sed -n ${i}p $in_map | awk '{print $1}'`
	read1=`sed -n ${i}p $in_map | awk '{print $2}'`
	read2=`sed -n ${i}p $in_map | awk '{print $3}'`
	out_bam_prefix=$tech"_to_"$ref_base
	header="@RG\tID:$tech\tLB:$tech\tPL:illumina\tSM:HB\tPU:$tech"
	echo "=== Start bwa align for $tech reads ==="
#:<<'END'

	# Align and merge
	echo "bwa mem -R$header -t$threads $ref $read1 | samtools view -hb -@$threads - > $out_bam_prefix.r1.bam"
	bwa mem -R$header -t$threads $ref $read1 | samtools view -hb -@$threads - > $out_bam_prefix.r1.bam
	echo ""
	
	echo "bwa mem -R$header -t$threads $ref $read2 | samtools view -hb -@$threads - > $out_bam_prefix.r2.bam"
	bwa mem -R$header -t$threads $ref $read2 | samtools view -hb -@$threads - > $out_bam_prefix.r2.bam
	echo ""

	# This merge does not really work...
	echo "samtools merge -@$threads -n -O bam $out_bam_prefix.bam $out_bam_prefix.r1.bam $out_bam_prefix.r2.bam"
	samtools merge -@$threads -n -O bam $out_bam_prefix.bam $out_bam_prefix.r1.bam $out_bam_prefix.r2.bam
	echo ""

	echo "samtools view -@$threads -q5 -hb $out_bam_prefix.bam > $out_bam_prefix.q5.bam"
	samtools view -@$threads -q5 -hb $out_bam_prefix.bam > $out_bam_prefix.q5.bam
	echo ""
#END

done
