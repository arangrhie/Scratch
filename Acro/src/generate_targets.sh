#!/bin/sh

if [[ $# -lt 2 ]]; then
  echo "Usage: $0 <annotate_regions.txt> <fasta>"
  echo "  annotate_regions.txt  a file with three columns: <name of the target region> <tab> <chr:start-end> <tab> <color in final annotation>"
  echo "  fasta                 reference fasta file to extract target regions. e.g. chm13v2.0.fa"
  exit 1
fi

target_regions=$1
fa=$2

module load samtools/1.21

LEN=`wc -l $target_regions | awk '{print $1}'`
for i in $(seq 1 $LEN)
do
  line=`sed -n ${i}p $target_regions`
  qry_name=`echo $line | awk '{print $1}'`
  region=`echo $line | awk '{print $2}'`
  if [[ $region == "." ]]; then
    continue
  fi

  if [[ -s $qry_name.fa ]]; then
    echo "$qry_name.fa exists. Skipping."
    continue
  fi
  
  direction=`echo $line | awk '{print $3}'`
  echo "Extracting $region from $fa to $qry_name.fa"

  if [[ $direction == "-" ]]; then
    opt="-i"
  else
    opt=""
  fi
  echo ">$qry_name" > $qry_name.fa
  samtools faidx $opt $fa $region | awk 'NR>1' >> $qry_name.fa
  samtools faidx $qry_name.fa
  
done
