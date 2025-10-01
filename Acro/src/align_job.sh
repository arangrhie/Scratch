#!/bin/bash

if [[ $# -lt 3 ]]; then
  echo "Usage: ./annotate_job.sh in_fa_list_prefix qry_fa len"
  echo "Check inputs. Exit."
  exit 0;
fi

list=$1.${SLURM_ARRAY_TASK_ID}
qry=$2
len=$3

for ref in $(cat $list)
do
  $tools/T2T-chrY/src/annotate.sh annotate_regions.txt $ref $qry
done


