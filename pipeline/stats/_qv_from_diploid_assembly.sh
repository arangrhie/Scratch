WEIGHT=0.5

PREFIX=10x_spnv_arrow


#awk '$1~/^10x_hap1/' $PREFIX.bayes.vcf > $PREFIX.bayes.hap1.vcf
#awk '$1~/^10x_hap2/' $PREFIX.bayes.vcf > $PREFIX.bayes.hap2.vcf

for HAP in hap1 hap2
do
NUM_SNP=`cat $PREFIX.bayes.$HAP.vcf |grep -v "#" | awk -F "\t" '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$8}' | tr ';' ' ' | sed s/AB=//g | awk -v WEIGHT=$WEIGHT '{if ($6 > WEIGHT) print $0}' | awk -v SUM=0 '{if (length($4) == length($5)) { SUM+=length($4); } else if (length($4) < length($5)) { SUM+=length($5)-length($4); } else { SUM+=length($4)-length($5)}} END { print SUM}'`

NUM_BP=`cat ${PREFIX}_hap1.numbp`
QV=`echo "$NUM_SNP $NUM_BP" |  awk '{print (-10*log($1/$2)/log(10))}'`

echo $NUM_SNP > ${PREFIX}_${HAP}_all.numsnp
echo $QV > ${PREFIX}_${HAP}_all.qv
done
