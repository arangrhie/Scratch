#!/bin/bash

syst=`uname -s`
arch=`uname -m`
name=`uname -n`

export BWA=/data/projects/phillippy/software/bwa/bwa-0.7.15
export PATH=$BWA:$PATH

if [ "$arch" = "x86_64" ] ; then
  arch="amd64"
fi
if [ "$arch" = "Power Macintosh" ] ; then
  arch="ppc"
fi

jobid=$SGE_TASK_ID
refid=`basename $1`
prefix=${refid/.fasta/}
map=${prefix}_map

if [ x$jobid = x -o x$jobid = xundefined -o x$jobid = x0 ]; then
   jobid=$2
fi

if test x$jobid = x; then
  echo Error: I need SGE_TASK_ID set, or a job index on the command line
  exit 1
fi

fofn=input/mapping.fofn

NUM_JOBS=`cat $fofn |wc -l |awk '{print $NF}'`
if [ $jobid -gt $NUM_JOBS ]; then
   echo "Invalid job id $jobid, max is $NUM_JOBS"
   exit
fi

reads=`cat $fofn |head -n $jobid |tail -n 1 |awk '{print $1}'`
readsPrefix=`basename $reads |sed s/.raw.fastq.gz//g |sed s/.fastq.gz//g |sed s/.fastq//g |sed s/.fq.gz//g |sed s/.fq//g`
prefix="${prefix}_${readsPrefix}"
type=`cat $fofn |head -n $jobid |tail -n 1 |awk '{print $NF}'`
echo "Running and will map reads $reads to $refid will prefix to $prefix.bam"

INPUT=""
OPTIONS=""
if [ $type == "pacbio" ]; then
  INPUT="$reads"
  OPTIONS="-x pacbio "
elif [ $type == "illumina" ]; then
   INPUT="${reads}_1.fastq.gz ${reads}_2.fastq.gz"
elif [ $type == "illuminaR1" ]; then
   readsV2=`echo $reads |sed s/_R1_/_R2_/g`
   INPUT="${reads} ${readsV2}"
else
   echo "Unknown data type $type"
   exit
fi

mkdir -p $map

if [ ! -e $map/$prefix.splitters.bam ]; then
#   if [ ! -e $map/$prefix.unsorted.bam ]; then
      echo "Runing bwa mem -t 16 $OPTIONS index/$refid $INPUT 2> logs/$prefix.err |samtools view -S -hb - > $map/$prefix.unsorted.bam"
      bwa mem -t 16 $OPTIONS index/$refid $INPUT 2> logs/$prefix.err |samtools view -Shb - > $map/$prefix.unsorted.bam
#   fi

   if [ ! -e $map/$prefix.splitters.unsorted.bam ]; then 
      echo "Extracting splitters/discordants"
      # extract
      samtools view -@16 -bF 1294 $map/$prefix.unsorted.bam > $map/$prefix.discordants.unsorted.bam
      samtools view -@16 -h $map/$prefix.unsorted.bam \
          | /data/projects/phillippy/software/lumpy-sv/scripts/extractSplitReads_BwaMem -i stdin \
          | samtools view -@16 -Shb - \
          > $map/$prefix.splitters.unsorted.bam
    fi
    if [ ! -e $map/$prefix.splitters.bam ]; then
      # sort
      echo "Sorting"
      ulimit -Su 160000
      ulimit -all

      samtools sort -@ 16 -m 2G -o $map/$prefix.bam -T $map/$prefix.tmp $map/$prefix.unsorted.bam
      samtools sort -@ 16 -m 2G -o $map/$prefix.discordants.bam -T $map/$prefix.tmp $map/$prefix.discordants.unsorted.bam
      samtools sort -@ 16 -m 2G -o $map/$prefix.splitters.bam -T $map/$prefix.tmp $map/$prefix.splitters.unsorted.bam
      samtools index $map/$prefix.bam
      samtools index $map/$prefix.discordants.bam
      samtools index $map/$prefix.splitters.bam
   fi
fi
