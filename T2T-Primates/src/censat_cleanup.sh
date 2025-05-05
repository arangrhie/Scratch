#!/bin/sh

# This script is used to clean up the censat 

IN=../input/all.chrom_breakdown.tsv
OUT=../input/all.chrom_breakdown.clean.tsv

echo "Sample_Chromosome\tCenSat\tSize" > $OUT

cat $IN | \
    awk -F"\t" 'BEGIN{OFS="\t"} \
	{if (NR>1 && $3>0) {print $1"_"$2, "HSat1A", $3}}' \
	>> $OUT

cat $IN | \
    awk -F"\t" 'BEGIN{OFS="\t"} \
	{if (NR>1 && $4>0) {print $1"_"$2, "HSat1B", $4}}' \
	>> $OUT

cat $IN | \
    awk -F"\t" 'BEGIN{OFS="\t"} \
	{if (NR>1 && $5>0) {print $1"_"$2, "HSat2", $5}}' \
	>> $OUT

cat $IN | \
    awk -F"\t" 'BEGIN{OFS="\t"} \
	{if (NR>1 && $6>0) {print $1"_"$2, "HSat3", $6}}' \
	>> $OUT

cat $IN | \
    awk -F"\t" 'BEGIN{OFS="\t"} \
	{if (NR>1 && $7>0) {print $1"_"$2, "bSat", $7}}' \
	>> $OUT

cat $IN | \
    awk -F"\t" 'BEGIN{OFS="\t"} \
	{if (NR>1 && $8>0) {print $1"_"$2, "Other", $8}}' \
	>> $OUT

cat $IN | \
    awk -F"\t" 'BEGIN{OFS="\t"} \
	{if (NR>1 && $9>0) {print $1"_"$2, "aSat", $9}}' \
	>> $OUT

cat $IN | \
    awk -F"\t" 'BEGIN{OFS="\t"} \
	{if (NR>1 && $10>0) {print $1"_"$2, "rDNA", $10}}' \
	>> $OUT

cat $IN | \
    awk -F"\t" 'BEGIN{OFS="\t"} \
	{if (NR>1 && $11>0) {print $1"_"$2, "ct", $11}}' \
	>> $OUT

cat $IN | \
    awk -F"\t" 'BEGIN{OFS="\t"} \
	{if (NR>1 && $12>0) {print $1"_"$2"_SD", "SD", $12}}' \
	>> $OUT

# Collect hap1
cat $OUT | \
	grep -v "hap2" | grep -v "PATERNAL" | grep -v "pat" \
	> ../input/all.chrom_breakdown.clean.hap1.tsv

# Collect hap2
cat $OUT | \
	grep -E "Sample|hap2|PATERNAL|pat" \
	> ../input/all.chrom_breakdown.clean.hap2.tsv

IN=../input/censat_summary.tsv
OUT=../input/censat_summary.clean.tsv

echo "Sample_Chromosome\tCenSat\tSize" > $OUT

cat $IN | \
    awk -F"\t" 'BEGIN{OFS="\t"} \
	{if (NR>1 && $2>0) {print $1, "HSat1A", $2}}' \
	>> $OUT

cat $IN | \
    awk -F"\t" 'BEGIN{OFS="\t"} \
	{if (NR>1 && $3>0) {print $1, "HSat1B", $3}}' \
	>> $OUT

cat $IN | \
    awk -F"\t" 'BEGIN{OFS="\t"} \
	{if (NR>1 && $4>0) {print $1, "HSat2", $4}}' \
	>> $OUT

cat $IN | \
    awk -F"\t" 'BEGIN{OFS="\t"} \
	{if (NR>1 && $5>0) {print $1, "HSat3", $5}}' \
	>> $OUT

cat $IN | \
    awk -F"\t" 'BEGIN{OFS="\t"} \
	{if (NR>1 && $6>0) {print $1, "bSat", $6}}' \
	>> $OUT

cat $IN | \
    awk -F"\t" 'BEGIN{OFS="\t"} \
	{if (NR>1 && $7>0) {print $1, "Other", $7}}' \
	>> $OUT

cat $IN | \
    awk -F"\t" 'BEGIN{OFS="\t"} \
	{if (NR>1 && $8>0) {print $1, "aSat", $8}}' \
	>> $OUT

cat $IN | \
    awk -F"\t" 'BEGIN{OFS="\t"} \
	{if (NR>1 && $9>0) {print $1, "rDNA", $9}}' \
	>> $OUT

cat $IN | \
    awk -F"\t" 'BEGIN{OFS="\t"} \
	{if (NR>1 && $10>0) {print $1, "ct", $10}}' \
	>> $OUT

cat $IN | \
    awk -F"\t" 'BEGIN{OFS="\t"} \
	{if (NR>1 && $11>0) {print $1"_SD", "SD", $11}}' \
	>> $OUT