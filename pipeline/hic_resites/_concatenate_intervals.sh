echo -e "FASTA\tRE\tInterval" > intervals.tds
cat *.interval >> intervals.tds 
awk '{print $1"_"$2"\t"$3}' intervals.tds > intervals.simple.tds
