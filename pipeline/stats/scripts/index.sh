#!/bin/sh

jobid=$SGE_TASK_ID
if [ x$jobid = x -o x$jobid = xundefined -o x$jobid = x0 ]; then
  jobid=$1
fi
if [ x$jobid = x ]; then
  echo Error: I need SGE_TASK_ID set, or a job index on the command line.
  exit 1
fi

file=`cat input/assembly.list | sort | head -n $jobid |tail -n 1`
asmName=`echo $file |awk -F "/" '{print $NF}'`
if [ ! -e index/$asmName.bwt ]; then
   echo "Need to index $file"
   mkdir -p index
   bwa index $file -p index/$asmName
fi
