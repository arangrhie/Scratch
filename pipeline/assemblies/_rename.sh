f1=$1
f2=$2

if [ -e $f1 ]; then
	mv $f1 $f2
	mv $f1.gaps $f2.gaps
	mv $f1.gaps.bed $f2.gaps.bed
	mv $f1.contigs.len $f2.contigs.len
	mv len/$f1.len len/$f2.len
fi
