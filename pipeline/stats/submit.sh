mkdir -p logs
#:<<'END'
for i in $(seq 102 105)  ; #    $(seq 11 11);
do
   file=`sed -n ${i}p input/assembly.list`
   prefix=`basename $file`
   prefix=${prefix/.fasta/}

   if [ ! -e $file.fai ]; then
      echo "qsub -cwd -V -j y -q phillippy.q -N index_"$prefix" -pe thread 1 -o logs/$prefix.index -S /bin/bash /home/rhiea/codes/_index.sh $file"
      qsub -cwd -V -j y -q phillippy.q -N index_"$prefix" -pe thread 1 -o logs/$prefix.index -S /bin/bash /home/rhiea/codes/_index.sh $file
      wait_for="-hold_jid index_"$prefix
   fi
#:<<'END'
   echo "qsub -cwd -V -j y -q phillippy.q $wait_for -N map_$prefix -t 1-11 -pe thread 16 -l mem_free=3g -o logs/$prefix.\\\$TASK_ID.out -S /bin/bash scripts/map.sh $file"
   qsub -cwd -V -j y -q phillippy.q $wait_for -N map_$prefix -t 1-11 -pe thread 16 -l mem_free=3g -o logs/$prefix.\$TASK_ID.out -S /bin/bash scripts/map.sh $file
#END
   echo "qsub -cwd -V -j y -q phillippy.q -N a_$prefix -hold_jid map_$prefix -pe thread 8 -l mem_free=20g,h_vmem=10G -o logs/$prefix.runAnalysis.out -S /bin/bash scripts/runAnalysis.sh $file"
   qsub -cwd -V -j y -q phillippy.q -N a_$prefix -hold_jid map_$prefix -pe thread 8 -l mem_free=10g -o logs/$prefix.runAnalysis.out -S /bin/bash scripts/runAnalysis.sh $file
   #qsub -cwd -V -j y -q phillippy.q,public.q -N a_$prefix -pe thread 8 -l mem_free=10g -o logs/$prefix.runAnalysis.out -S /bin/bash scripts/runAnalysis.sh $file

   echo "qsub -cwd -V -j y -q phillippy.q -N qv_$prefix -pe thread 16 -l mem_free=20g,h_vmem=20G -o logs/$prefix.runQV.out -S /bin/bash scripts/runQV.sh $file"
#END

#:<<'END'
   echo "Busco"
   echo "qsub -cwd -V -j y -q phillippy.q -N b_$prefix -pe thread 12 -o logs/$prefix.busco.out -S /bin/bash scripts/runBusco.sh ../$file 12"
   qsub -cwd -V -j y -q phillippy.q -N b_$prefix -pe thread 12 -o logs/$prefix.busco.out -S /bin/bash scripts/runBusco.sh ../$file 12
#END
done
#END

