threads=$1
mem_per_threads=2000000000

#echo "/data/projects/phillippy/software/jellyfish-2.2.6/bin/jellyfish count -C -m 21 -s $mem_per_threads -t $threads 10x/Hummingbird_S1_L001_R[1-2]_001.fastq -o k21_10xreads.jf"
#/data/projects/phillippy/software/jellyfish-2.2.6/bin/jellyfish count -C -m 21 -s $mem_per_threads -t $threads 10x/Hummingbird_S1_L001_R[1-2]_001.fastq -o k21_reads.jf

echo "/data/projects/phillippy/software/jellyfish-2.2.6/bin/jellyfish histo -t $threads k21_reads.jf > reads.histo"
/data/projects/phillippy/software/jellyfish-2.2.6/bin/jellyfish histo -t $threads k21_reads.jf > reads.histo

