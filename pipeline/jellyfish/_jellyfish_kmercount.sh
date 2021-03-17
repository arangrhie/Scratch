k=$1
threads=$2
fasta=$3

/data/projects/phillippy/software/jellyfish-2.2.6/bin/jellyfish count -m $k -s 1G -t $threads $fasta -o ${fasta/.fasta/.$k.jf}
