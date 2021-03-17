assembly=$1
size=$2

mkdir -p len
mkdir -p stats

if [ ! -e "len/$assembly.len" ]; then
	java -jar -Xmx4g /home/rhiea/codes/fastaContigSize.jar $assembly
	mv $assembly.len len/$assembly.len
fi

if [ "$size" = "" ]; then
	java -jar -Xmx4g /home/rhiea/codes/lenCalcNGStats.jar len/$assembly.len 1100000000 > stats/$assembly.len.stat
else
	java -jar -Xmx4g /home/rhiea/codes/lenCalcNGStats.jar len/$assembly.len $size > stats/$assembly.len.stat
fi
cat stats/$assembly.len.stat
echo "N-bases: "`awk '{sum+=$2; sumN+=$3} END {print (sum-sumN)}' len/$assembly.len`

./_split_scaffolds_by_gaps.sh $assembly

./_get_gaps_stat.sh $assembly

