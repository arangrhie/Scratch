## For calculating actual coverage of reads
cd aligned/
#shuf -n 2000000 merged_nodups.txt.subsample > merged_nodups.txt.subsample.2M	# for arima and dovetail chicago
shuf -n 2000000 merged_nodups.txt > merged_nodups.txt.subsample.2M		# for phase and dovetail hiC (2)
java -jar -Xmx1g /home/rhiea/codes/juicerToBed.jar merged_nodups.txt.subsample.2M > merged_nodups.txt.subsample.2M.bed
java -jar -Xmx1g /home/rhiea/codes/bedToDistance.jar merged_nodups.txt.subsample.2M.bed > merged_nodups.txt.subsample.2M.bed.dist
cd ../

bed=aligned/merged_nodups.txt.subsample.2M.bed
len=/data/projects/phillippy/projects/hummingbird/juicer_hic/pacbio_falcon_pa_pilon_ponly.fasta.len.bed
cov_hist=cov_hist_rand_2M
mkdir -p $cov_hist

echo "=== Get coverage by insert size ==="
for isize in 500 1000 5000 10000 50000 100000 500000 1000000 5000000 10000000 20000000
do
        echo "... > $isize"
        #awk -v isize=$isize '($3-$2) >= isize' $bed.sort | bedtools coverage -hist -a - -b $len > $cov_hist/$isize.bed
        awk -v isize=$isize '$NF >= isize' $bed | bedtools coverage -hist -a - -b $len > $cov_hist/$isize.bed
        awk -v isize=$isize '$1=="all" && $2>=1 {sum += $NF} END {print isize"\t"sum}' $cov_hist/$isize.bed >> cov_over_1.2M.bed
done

