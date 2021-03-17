in_prefix=$1
ref=$2

in_prefix=$1
ref=$2

if [ ! -e "$ref.fai" ] ; then
	echo "samtools dict -o ${ref/.fasta/.dict} $ref"
	samtools dict -o ${ref/.fasta/.dict} $ref
	echo "samtools faidx $ref"
	samtools faidx $ref
fi

# Need to provide unaligned.bam
java -jar -Xmx24g /data/projects/phillippy/software/picard-tools-2.7.1/picard.jar MergeBamAlignment \
	UNMAPPED=$in_prefix.unmapped.bam \
	R1_ALIGNED=$in_prefix.r1.bam \
	R2_ALIGNED=$in_prefix.r1.bam \
	R=$ref \
	O=$in_prefix.mrg.bam

