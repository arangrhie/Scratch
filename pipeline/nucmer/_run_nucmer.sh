
ref=$1
qry=$2

qsub -V -q phillippy.q -P large_mem -j y -cwd -N nucmer -pe thread 4 -l mem_free=60g -o logs/nucmer.${qry/.fasta/}_to_${ref/.fasta/}.log -S /bin/bash _nucmer.sh $ref $qry
