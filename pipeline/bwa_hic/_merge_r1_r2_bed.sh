prefix=$1
awk '{print $0"\tR1"}' $prefix.r1.bed > $prefix.R1.bed
awk '{print $0"\tR2"}' $prefix.r2.bed > $prefix.R2.bed
sort -k4 $prefix.R1.bed $prefix.R2.bed > $prefix.merged.bed

