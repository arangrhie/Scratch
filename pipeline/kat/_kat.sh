assembly=$1
threads=$2

out=`echo $assembly | awk -F "/" '{print $NF}'`
out=${out/.fasta/}
mkdir -p $out
out=$out/kat_comp
#options="-H 10000000000 -I 10000000000 -m 27 -h"
options="-H 10000000000 -I 10000000000 -m 22 -h"

echo "/data/projects/phillippy/software/kat-2.3.4/bin/kat comp -o $out -t $threads $options -O 'R1.BC.fastq R2.BC.fastq' $assembly"
/data/projects/phillippy/software/kat-2.3.4/bin/kat comp -o $out -t $threads $options -O 'R1.BC.fastq R2.BC.fastq' $assembly

