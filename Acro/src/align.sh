#!/bin/sh

if [[ $# -lt 4 ]]; then
  echo "Usage: align.sh ref.fa qry.fa len"
  echo "Uses minimap2 to map ref to qry, then inverse coordinates and filter for 10kb alignments."
  echo "  ref.fa    reference fasta file"
  echo "  qry.fa    querey fasta file"
  echo "  len       minimum length of alignment to keep"
  exit 0
fi

ref=$1 # assembly fasta
ref_name=`basename $ref`
ref_name=`echo $ref_name | sed 's/.gz$//g' | sed 's/\.fa$//g' | sed 's/\.fasta$//g'`

qry=$2 # region of interest fasta
qry_name=`basename $qry`
qry_name=`echo $qry_name | sed 's/.gz$//g' | sed 's/\.fa$//g' | sed 's/\.fasta$//g'`

len=$3

if [[ -s ${qry_name}_to_${ref_name}.filter.bed ]]; then
  echo "${qry_name}_to_${ref_name}.filter.bed already exists. Exit."
  exit 0
fi

echo $qry_name

module load samtools/1.21
module load minimap2/2.28
module load rustybam/0.1.33

set -x
minimap2 \
    -x asm10 \
    --eqx \
    --MD \
    -t ${SLURM_CPUS_PER_TASK} \
    -c \
    $qry \
    $ref \
    > ${ref_name}_to_${qry_name}.paf

rustybam invert ${ref_name}_to_${qry_name}.paf > ${qry_name}_to_${ref_name}.paf
set +x

awk -v len=$len '{ if( $11 >= len ) print }' ${qry_name}_to_${ref_name}.paf > ${qry_name}_to_${ref_name}.filter.paf
cat ${qry_name}_to_${ref_name}.filter.paf |\
  awk -v OFS='\t' -v qry_name=${qry_name} '{if ($10<$11) {idy=100*$10/$11;} else {idy=100*$11/$10;} print $6,$8,$9,qry_name,idy,$5}' |\
  sort -k1,1V -k2,2n - > ${qry_name}_to_${ref_name}.filter.bed

echo "Final alignment file: ${qry_name}_to_${ref_name}.filter.bed, filtered at > $len bp"


