#if [ ! -e $1.contigs.len ]; then
java -jar -Xmx4g $codes/fastaGetGaps.jar $1 $1.gaps
awk -F "\t" '$4>3 {print $1"\t"$2"\t"$3}' $1.gaps > $1.gaps.bed
awk '{print $1"\t0\t"$(NF-1)}' len/$1.len > len/$1.len.bed
bedtools subtract -a len/$1.len.bed -b $1.gaps.bed | awk '{print $NF-$(NF-1)}' > $1.contigs.len
#fi
java -jar -Xmx1g /home/rhiea/codes/lenCalcNGStats.jar $1.contigs.len 1100000000 1
