ref=$1	# ref.fasta
qry=$2	# qry.fasta

ref_prefix=`basename $ref`
ref_prefix=${ref_prefix/.fasta/}
qry_prefix=`basename $qry`
qry_prefix=${qry_prefix/.fasta/}

prefix=${qry_prefix}_to_${ref_prefix}
mkdir -p $prefix
out=$prefix/$prefix

#:<<'END'
echo "/data/projects/phillippy/software/MUMmer3.23/nucmer -p=$out -l 100 -c 1000 $ref $qry"
/data/projects/phillippy/software/MUMmer3.23/nucmer -p=$out -l 100 -c 1000 $ref $qry

echo "/data/projects/phillippy/software/MUMmer3.23/dnadiff -p $out.dnadiff -d $out.delta"
/data/projects/phillippy/software/MUMmer3.23/dnadiff -p $out.dnadiff -d $out.delta
#END

echo "/data/projects/phillippy/software/MUMmer3.23/mummerplot --large --fat -t png -p $out.dnadiff.1delta.plot $out.dnadiff.1delta"
/data/projects/phillippy/software/MUMmer3.23/mummerplot --large --fat -t png -p $out.dnadiff.1delta.plot $out.dnadiff.1delta

