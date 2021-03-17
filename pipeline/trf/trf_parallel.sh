
jobid=$SGE_TASK_ID

if test x$jobid = x; then
  echo Error: I need SGE_TASK_ID set, or a job index on the command line
  exit 1
fi

PREFIX=$1
LEN=$PREFIX.lens

line=`cat $LEN |head -n $jobid |tail -n 1`
contig=`echo $line |awk '{print $1}'`

echo "/home/rhiea/tools/trf409.legacylinux64 split_fasta/$contig.fa 2 7 7 80 10 50 500 -f -d -m -h"
/home/rhiea/tools/trf409.legacylinux64 split_fasta/$contig.fa 2 7 7 80 10 50 500 -f -d -m -h

