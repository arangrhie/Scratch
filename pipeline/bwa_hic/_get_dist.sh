bam=$1

bedtools bamtobed -i $bam | head -n 100000 > ${bam/.bam/.bed}
java -jar -Xmx1g /home/rhiea/codes/bedToDistance.jar ${bam/.bam/.bed} > ${bam/.bam/.dist}
head -n 50000 ${bam/.bam/.dist} > ${bam/.bam/.dist}.5000
