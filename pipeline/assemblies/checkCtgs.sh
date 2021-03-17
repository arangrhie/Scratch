for file in `ls *.fasta`; do
   echo "***** $file ****"
   cat $file |grep ">"|awk '{print $1}' |sort |uniq -c |awk '{print $1}'|sort |uniq -c
   echo "****************"
done
