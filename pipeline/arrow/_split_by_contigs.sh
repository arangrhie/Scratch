prefix=$1
if [ ! -e $prefix.fasta.len ] ; then
	echo "No $prefix.fasta.len found. Create one..."
	echo "java -jar -Xmx1g /home/rhiea/codes/fastaContigSize.jar $prefix.fasta"
	java -jar -Xmx1g /home/rhiea/codes/fastaContigSize.jar $prefix.fasta
fi
NUM_JOBS=`wc -l $prefix.fasta.len | awk '{print $1}'`

echo "dataset create --type AlignmentSet $prefix/$prefix.xml $prefix/[0-9]*.aln.bam"
dataset create --type AlignmentSet $prefix/$prefix.xml $prefix/[0-9]*.aln.bam

mkdir -p $prefix/by_contigs
echo "dataset split --contig --maxChunks $NUM_JOBS --chunks $NUM_JOBS --outdir $prefix/by_contigs $prefix/$prefix.xml"
dataset split --contig --maxChunks $NUM_JOBS --chunks $NUM_JOBS --outdir $prefix/by_contigs $prefix/$prefix.xml

