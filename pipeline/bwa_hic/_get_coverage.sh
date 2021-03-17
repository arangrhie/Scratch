bam_prefix=$1	# expect to have <$bam_prefix>.r1.bam and >$bam_prefix>.r2.bam
ref_fa=$2	# full path to reference fasta

#:<<END
if [ ! -e "$bam_prefix.r1.bed" ] ; then
echo "bedtools bamtobed -split -i $bam_prefix.r1.bam > $bam_prefix.r1.bed"
bedtools bamtobed -split -i $bam_prefix.r1.bam > $bam_prefix.r1.bed
echo ""
echo "bedtools bamtobed -split -i $bam_prefix.r2.bam > $bam_prefix.r2.bed"
bedtools bamtobed -split -i $bam_prefix.r2.bam > $bam_prefix.r2.bed
echo ""
fi
#END

if [ ! -e "$ref_fa.len.bed" ] ; then
	if [ ! -e "$ref_fa.len" ] ; then
		echo "java -jar -Xmx1g /home/rhiea/codes/fastaContigSize.jar $ref_fa"
		java -jar -Xmx1g /home/rhiea/codes/fastaContigSize.jar $ref_fa
	fi
	echo "awk '{print $1"\t0\t"$NF}' $ref_fa.len > $ref_fa.len.bed"
	awk '{print $1"\t0\t"$NF}' $ref_fa.len > $ref_fa.len.bed
fi

if [ ! -e "$bam_prefix.r1.cov" ] ; then
echo "bedtools coverage -a $bam_prefix.r1.bed -b $ref_fa.len.bed > $bam_prefix.r1.cov"
bedtools coverage -a $bam_prefix.r1.bed -b $ref_fa.len.bed > $bam_prefix.r1.cov
echo ""

echo "bedtools coverage -a $bam_prefix.r2.bed -b $ref_fa.len.bed > $bam_prefix.r2.cov"
bedtools coverage -a $bam_prefix.r2.bed -b $ref_fa.len.bed > $bam_prefix.r2.cov
echo ""
fi
#END

if [ ! -e "$bam_prefix.r1.count_reads" ] ; then
echo "java -jar -Xmx4g /home/rhiea/codes/bedCountSplits.jar $bam_prefix.r1.bed $bam_prefix.r1.count_reads"
java -jar -Xmx4g /home/rhiea/codes/bedCountSplits.jar $bam_prefix.r1.bed $bam_prefix.r1.count_reads

echo "java -jar -Xmx4g /home/rhiea/codes/bedCountSplits.jar $bam_prefix.r2.bed $bam_prefix.r2.count_reads"
java -jar -Xmx4g /home/rhiea/codes/bedCountSplits.jar $bam_prefix.r2.bed $bam_prefix.r2.count_reads
fi

