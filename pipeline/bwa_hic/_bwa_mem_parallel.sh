ref=$1
ref_base=${ref/.fasta/}
ref_base=${ref/.fa/}
ref_base=${ref/.fna/}

i=$2
threads=$3

in_map=hic_reads.map

export BWA=/data/projects/phillippy/software/bwa/bwa-0.7.15
export PATH=$BWA:$PATH


tech=`sed -n ${i}p $in_map | awk '{print $1}'`
read1=`sed -n ${i}p $in_map | awk '{print $2}'`
read2=`sed -n ${i}p $in_map | awk '{print $3}'`
out_bam_prefix=$tech"_to_"$ref_base.paired
header="@RG\tID:$tech\tLB:$tech\tPL:illumina\tSM:HB\tPU:$tech"
echo "=== Start bwa align for $tech reads ==="

# Align and merge
echo "bwa mem -R$header -t$threads $ref $read1 $read2 | samtools view -hb -@$threads - > $out_bam_prefix.bam"
bwa mem -R$header -t$threads $ref $read1 $read2 | samtools view -hb -@$threads - > $out_bam_prefix.bam
echo ""

echo "samtools view -@$threads -q5 -hb $out_bam_prefix.bam > $out_bam_prefix.q5.bam"
samtools view -@$threads -q5 -hb $out_bam_prefix.bam > $out_bam_prefix.q5.bam
echo ""

echo "samtools stats -i 1000000000 $out_bam_prefix.bam > $out_bam_prefix.stats"
samtools stats -i 1000000000 $out_bam_prefix.bam > $out_bam_prefix.stats
echo ""
