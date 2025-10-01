# Final issues in v1.1

Merge issues
```sh
module load bedtools
ver=v1.1
hifi_issues=/data/CARD_AUX/KOLF/polish/mapping/v1.1.dip_hifi/v1.1_dip.hifi.pri.issues.bed
ont_issues=/data/CARD_AUX/KOLF/polish/mapping/v1.1.dip_ont_ul/v1.1_dip.ont.pri.issues.bed

gap=../../pattern/v1.1.exclude.bed

cat ../../rDNA/45S_to_v1.1.mashmap.bed | awk '$(NF-1)>90' | bedtools merge -s -d 100000 -i - > rDNA.bed
rDNA=rDNA.bed

# extend rDNA.bed to flank gap
cut -f1-3 $gap $rDNA | sort -k1,1V -k2,2n | bedtools merge -d 100000 -i - | bedtools subtract -a - -b $gap > rDNA_refined.bed


# chromosome ends
sizes=../pattern/$ver.fa.fai
grep -v "chrM" $sizes | awk '{print $1"\t0\t100\n"$1"\t"($2-100)"\t"$2}' > ${ver}.end.bed

bedtools intersect -u -a $hifi_issues -b $ont_issues | grep -v "chrM" | grep -v "chrY" > ${ver}.issues.bed
bedtools merge -d 100000 -i ${ver}.issues.bed > ${ver}.issues.100k_mrg.bed
cat ${ver}.end.bed rDNA_refined.bed $gap | cut -f1-3 | sort -k1,1 -k2,2n | bedtools intersect -u -a ${ver}.issues.100k_mrg.bed -b - > ${ver}.issues.100k_mrg.rDNA_end_gap.bed
bedtools subtract -A -a ${ver}.issues.bed -b ${ver}.issues.100k_mrg.rDNA_end_gap.bed > ${ver}.issues.no_rDNA.no_end.no_gap.bed

bed=$gap
gap=`awk '{sum+=$3-$2} END {print sum}' $bed`

bed=rDNA_refined.bed
rDNA=`awk '{sum+=$3-$2} END {print sum}' $bed`

bed=${ver}.issues.no_rDNA.no_end.no_gap.bed
assembly=`awk '{sum+=$3-$2} END {print sum}' $bed`


```
