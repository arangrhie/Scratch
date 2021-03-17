# Rename haplotype fasta if naming conflicts
hap1=$1
hap1=${hap1/.fasta/}
hap2=$2
hap2=${hap2/.fasta/}

merged=$3
merged=${merged/.fasta}

:<<'END'
java -jar -Xmx1g /home/rhiea/codes/fastaContigSize.jar $hap1.fasta
cut -f1 $hap1.fasta.len > $hap1.fasta.list
awk '{print $0"_hap1"}' $hap1.fasta.list > $hap1.fasta.newname.list
java -jar -Xmx1g /home/rhiea/codes/fastaExtractFromList.jar $hap1.fasta $hap1.fasta.list $hap1.newname.fasta $hap1.fasta.newname.list

java -jar -Xmx1g /home/rhiea/codes/fastaContigSize.jar $hap2.fasta
cut -f1 $hap2.fasta.len > $hap2.fasta.list
awk '{print $0"_hap2"}' $hap2.fasta.list > $hap2.fasta.newname.list
java -jar -Xmx1g /home/rhiea/codes/fastaExtractFromList.jar $hap2.fasta $hap2.fasta.list $hap2.newname.fasta $hap2.fasta.newname.list

cat $hap1.newname.fasta $hap2.newname.fasta > $merged.fasta
END

rm $hap1.fasta.*list $hap2.fasta.*list
mv $hap1.newname.fasta $hap1.fasta
mv $hap2.newname.fasta $hap2.fasta
