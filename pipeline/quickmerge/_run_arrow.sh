#!/bin/bash
## _run_pbalign.sh

reference=$1	# .fasta assembly
prefix=${reference/.fasta/}
BAM_LIST="subreads.bam.list"
LEN=`wc -l $BAM_LIST | awk '{print $1}'`


if [ ! -e $prefix.fasta.len ] ; then
	echo "$prefix.fasta.len not found. Generating len file..."
	java -jar -Xmx1g /home/rhiea/codes/fastaContigSize.jar $prefix.fasta
fi

NUM_CONTIGS=`wc -l $prefix.fasta.len | awk '{print $1}'`

:<<'END'
mkdir -p logs
qsub -V -q phillippy.q,public.q -j y -cwd -N pbalign.$prefix -tc 60 -t 1-$LEN -pe thread 8 -l mem_free=6g -o logs/pbalign.$prefix.\$TASK_ID.log -S /bin/bash _pbalign_parallel.sh $reference $BAM_LIST
#qsub -V -q phillippy.q,public.q -j y -cwd -N pbalign.$prefix -tc 20 -t 1-$LEN -pe thread 8 -l mem_free=6g -o logs/pbalign.$prefix.\$TASK_ID.log -S /bin/bash _pbalign_parallel.sh $reference $BAM_LIST
END

#qsub -V -q phillippy.q,public.q -j y -hold_jid pbalign.$prefix -cwd -N split.$prefix -pe thread 2 -o logs/split.$prefix.log -S /bin/bash _split_by_contigs.sh $prefix
#qsub -V -q phillippy.q,public.q -j y -cwd -N split.$prefix -pe thread 2 -o logs/split.$prefix.log -S /bin/bash _split_by_contigs.sh $prefix

if [ ! -e $reference.fai ] ; then
	echo "$reference.fai not found. Start indexing..."
	samtools faidx $reference
fi
#END
#echo "qsub -V -q phillippy.q,public.q -j y -cwd -hold_jid split.$prefix -N arrow.$prefix -tc 40 -t 1-$NUM_CONTIGS -pe thread 24 -l mem_free=6g -o logs/arrow.$prefix.\$TASK_ID.log -S /bin/bash _consensus_parallel.sh $prefix"
qsub -V -q phillippy.q,public.q -j y -cwd -N arrow.$prefix -tc 40 -t 1-$NUM_CONTIGS -pe thread 12 -l mem_free=6g -o logs/arrow.$prefix.\$TASK_ID.log -S /bin/bash _consensus_parallel.sh $prefix
#qsub -V -q phillippy.q,public.q -j y -cwd -hold_jid split.$prefix -N arrow.$prefix -tc 40 -t 1-$NUM_CONTIGS -pe thread 12 -l mem_free=6g -o logs/arrow.$prefix.\$TASK_ID.log -S /bin/bash _consensus_parallel.sh $prefix
#qsub -V -q phillippy.q -P large_mem -j y -cwd -N arrow -tc 40 -t 8-$NUM_CONTIGS -pe thread 8 -o logs/arrow.$prefix.\$TASK_ID.log -S /bin/bash _consensus_parallel.sh $prefix

qsub -V -q phillippy.q,public.q -j y -cwd -hold_jid arrow.$prefix -N merge.$prefix -pe thread 1 -o logs/arrow.$prefix.\$TASK_ID.log -S /bin/bash _merge.sh $prefix

