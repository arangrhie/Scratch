fasta=$1	# msNA12878_pseudohap2.1.fasta

name=${fasta/.1.fasta/}
if ! [ -e $fasta.len ]; then
	java -jar -Xmx1g /home/rhiea/codes/fastaContigSize.jar $fasta
fi
awk '{print $1}' $fasta.len > $name.name
awk '{print $1"_hap1"}' $fasta.len > $name.name.hap1
awk '{print $1"_hap2"}' $fasta.len > $name.name.hap2
java -jar -Xmx2g /home/rhiea/codes/fastaExtractFromList.jar $name.1.fasta $name.name ${name}_hap1.fasta $name.name.hap1
java -jar -Xmx2g /home/rhiea/codes/fastaExtractFromList.jar $name.2.fasta $name.name ${name}_hap2.fasta $name.name.hap2

