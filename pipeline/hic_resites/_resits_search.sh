
assembly=$1
assembly=${assembly/.fasta/}
assembly=${assembly/.fa/}
assembly=${assembly/.fna/}

resite="GANTC"
if [ ! -e "$assembly.$resite.bed" ] ; then
	java -jar -Xmx4g /home/rhiea/codes/fastaFindSequence.jar $resite $assembly.fasta > $assembly.$resite.bed
fi
:<<'END'
if [ ! -e "$assembly.len" ]; then
        java -jar -Xmx1g /home/rhiea/codes/fastaContigSize.jar $1
        awk '{print $1"\t0\t"$3}' $1.len > $assembly.len
fi
if [ ! -e "$assembly.$resite.contigs.len" ] ; then
	cut -f 1 $assembly.$resite.bed | sort -u > $assembly.$resite.contig.list
	java -jar -Xmx1g /home/rhiea/codes/txtContains.jar $assembly.len $assembly.$resite.contig.list 1 > $assembly.$resite.contigs.len
fi
if [ ! -e "$assembly.$resite.v.bed" ] ; then
	bedtools subtract -a $assembly.$resite.contigs.len -b $assembly.$resite.bed > $assembly.$resite.v.bed
	awk -v assembly=$assembly -v resite=$resite '{print assembly"\t"resite"\t"$3-$2}' $assembly.$resite.v.bed > $assembly.$resite.interval
fi
END
echo "Finished to generate $assembly.$resite.bed"

# echo "Finished to generate $assembly.$resite.interval"

