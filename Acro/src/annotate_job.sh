#!/bin/bash

if [[ $# -lt 2 ]]; then
  echo "Usage: ./annotate_job.sh chr_assigned_fa.map out chr"
  echo "  chr_assigned_fa.map  prefix of the <sample> and <path/to/input/fasta> mapping file"
  echo "                       e.g., chr_assigned_fa.map for chr_assigned_fa.map.10, chr_assigned_fa.map.11, etc."
  echo "  out                  output directory"
  echo "  chr                  reference chromosome to fix the annotation to. optional."
  echo "Check inputs. Exit."
  exit 0;
fi

src=~/Scratch/Acro # change this to your own path to Scratch/Acro

map=$1.${SLURM_ARRAY_TASK_ID}
if [[ ! -s $map ]]; then
  echo "Using original map $1"
  map=$1
fi

out=$2
chr=$3

LEN=$(wc -l $map | awk '{print $1}')

for i in $(seq 1 $LEN)
do
  line=$(sed -n "${i}p" $map)
  sample=$(echo $line | awk '{print $1}')
  asm=$(echo $line | awk '{print $2}')
  $src/src/annotate.sh $src/input/annotate_regions.txt $sample $asm $out $chr
done


