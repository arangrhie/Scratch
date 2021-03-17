for file in `ls *.fasta`; do
   for file2 in `ls *.fasta`; do
      fID=`ls *.fasta |grep -n "^$file$" |awk -F ":" '{print $1}'`
      f2ID=`ls *.fasta |grep -n "^$file2$" | awk -F ":" '{print $1}'`

      if [ $fID -lt $f2ID ]; then
         if [ $file != $file2 ]; then
            prefix=`echo $file $file2 |sed s/.fasta//g |awk '{print $1"_"$2}'`
            if [ ! -e $prefix.report ]; then
               echo "/data/projects/phillippy/software/MUMmer3.23/nucmer -mumref -l 100 -c 1000 -p $prefix $file $file2"
               echo "/data/projects/phillippy/software/MUMmer3.23/dnadiff -p $prefix -d $prefix.delta"
               echo "/data/projects/phillippy/software/MUMmer3.23/mummerplot --large --fat --png $prefix.1delta"
            fi
         fi
      fi
   done
done
