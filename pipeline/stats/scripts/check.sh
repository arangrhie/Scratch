for file in `ls map/*.bam`; do
   echo $file
   samtools view $file |head > /dev/null
done
