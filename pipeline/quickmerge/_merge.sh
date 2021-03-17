prefix=$1
cat $prefix/out/*.fasta > ${prefix}_arrow.fasta
cat $prefix/out/*.fastq > ${prefix}_arrow.fastq
cat $prefix/out/*.gff > ${prefix}_arrow.gff
sed -i 's/|arrow/_arrow/g' ${prefix}_arrow.fasta
sed -i 's/|arrow/_arrow/g' ${prefix}_arrow.fastq
sed -i 's/|arrow/_arrow/g' ${prefix}_arrow.gff

echo "Check the ${prefix}_arrow.fasta file, then remove $prefix directory"
