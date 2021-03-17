fasta=$1	# <any>.fasta
new_prefix=$2
name=${fasta/.fasta/}

if ! [ -e $fasta.len ]; then
	java -jar -Xmx1g /home/rhiea/codes/fastaContigSize.jar $fasta
fi
awk '{print $1}' $fasta.len > $name.name
awk -v prefix=$new_prefix '{print prefix"_"$1}' $fasta.len | sed 's/|/_/g' > $name.name.new
java -jar -Xmx2g /home/rhiea/codes/fastaExtractFromList.jar $fasta $name.name $new_prefix"_"$fasta $name.name.new

