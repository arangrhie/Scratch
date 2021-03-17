cat aligned/merged_nodups.txt |awk '{if ($2 == $6 && $1 ==0 && $5 == 0) print $0}'|awk '{if ($3 < $7) print $7-$3; else print $3-$7}' > aligned/merged_nodups.dist
