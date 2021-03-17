assembly_fa=$1	#.fasta
prefix=${assembly_fa/.fasta/}
pilon_dir=${prefix}_pilon

mkdir -p $pilon_dir

echo "=== Create symbolic links ==="
ln -s /data/projects/phillippy/software/PilonGrid/pilon_10X.sh $pilon_dir
ln -s $PWD/$prefix/pos_sorted_bam.bam $pilon_dir/pos_sorted_bam.bam
ln -s $PWD/$prefix/pos_sorted_bam.bam.bai $pilon_dir/pos_sorted_bam.bai
ln -s $PWD/$assembly_fa $pilon_dir
echo ""

echo "cd $pilon_dir"
cd $pilon_dir
echo "ls"
ls
echo ""

echo "=== Start running pilon_10X.sh ==="
echo "./pilon_10X.sh $assembly_fa"
./pilon_10X.sh $assembly_fa

