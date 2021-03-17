nucmer -mumref -l 100 -c 1000 -banded -D 10
delta-filter -i 95 -o 95 out.delta > out.best.delta
dnadiff -d out.best.delta
