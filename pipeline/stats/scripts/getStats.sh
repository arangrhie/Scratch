> errors.csv
FIRST=1

for file in `ls *_Features.txt`; do
   name=`echo $file |sed s/_Features.txt//g`
   vcf="$name.lumpy.vcf"

   if [ $FIRST -eq 1 ]; then
      list=`cat *_Features.txt  |awk '{print $2}'|sort |uniq -c |awk '{print $NF}' |tr '\n' ','`
      lumpyList=`cat *.lumpy.vcf  |grep -v "#"  |awk -F "\t" '{print $8}'|awk -F ";" '{print $1}'|sort |uniq -c |awk '{print $NF}' |tr '\n' ','`
      FIRST=0
      echo "ASM,${list}${lumpyList}" > errors.csv
   fi

   l=""
   for feature in `echo $list |tr ',' ' '`; do
      t=`grep -c $feature $file`
      if [ x$l == "x" ]; then
         l=$t
      else
         l="$l,$t"
     fi
   done
   for feature in `echo $lumpyList |tr ',' ' '`; do
      t=`cat $vcf 2>/dev/null |grep -v "#"  |awk -F "\t" '{print $8}'|awk -F ";" '{print $1}'|sort |uniq -c |grep $feature |awk '{print $1}'`
      l="$l,$t"
   done
   echo "$name,${l}" >> errors.csv
done
tail -n 6 *.busco/run_SAMPLE/short_summary_SAMPLE |grep -v "^$" |awk '{if (match($0, "busco")) { NAME=substr($2, 1, index($2, "/")-1); COUNT=0; } else if (COUNT < 5) { NAME=NAME"\t"$1; COUNT++;} else if (COUNT==5) {print NAME" "$1}}' > busco.csv
tail *.qv|grep -v "^$" |awk '{if (match($0, "qv")) { NAME=$2;  } else { print NAME"\t"$0}}' > qv.csv

java GetFastaStats -o -min 0 -genomeSize 1100000000 *.fasta > scf.stats 2>&1
java GetFastaStats -o -min 0 -split NNN -genomeSize 1100000000 *.fasta > ctg.stats 2>&1

cat scf.stats |awk '{if (match($1, "Processing")) { NAME=$NF; } else if (match($0, "Total units")) { NUM=$NF; } else if (match($1, "BasesInFasta")) { SUM=$NF; } else if (match($1, "Max")) { MAX=$NF; } else if (match($1, "NG25")) { NG25=$2; } else if (match($1, "NG50")) { NG50=$2; } else if (match($1, "NG75")) { print NAME"\t"SUM"\t"NUM"\t"MAX"\t"NG25"\t"NG50"\t"$2}}' > scf.csv
cat ctg.stats |awk '{if (match($1, "Processing")) { NAME=$NF; } else if (match($0, "Total units")) { NUM=$NF; } else if (match($1, "BasesInFasta")) { SUM=$NF; } else if (match($1, "Max")) { MAX=$NF; } else if (match($1, "NG25")) { NG25=$2; } else if (match($1, "NG50")) { NG50=$2; } else if (match($1, "NG75")) { print NAME"\t"SUM"\t"NUM"\t"MAX"\t"NG25"\t"NG50"\t"$2}}' > ctg.csv

