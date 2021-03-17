## _pbalign_parallel.sh

reference=$1	# .fasta assembly
prefix=${reference/.fasta/}

jobid=$SGE_TASK_ID
if test x$jobid = x ; then
  echo Error: I need SGE_TASK_ID set, or a job index on the command line
  exit 1
fi

BAM_LIST=$2
line=`sed -n ${jobid}p $BAM_LIST`

echo "Mapping $line to $reference"
mkdir -p tmpdir
mkdir -p $prefix

echo "pbalign --tmpDir=`pwd`/tmpdir --minAccuracy=0.75 --minLength=50 --minAnchorSize=12 --maxDivergence=30 --concordant --algorithm=blasr --algorithmOptions=--useQuality --maxHits=1 --hitPolicy=random --seed=1 --nproc=8 $line $reference $prefx/$jobid.aln.bam"
pbalign --tmpDir=`pwd`/tmpdir --minAccuracy=0.75 --minLength=50 --minAnchorSize=12 --maxDivergence=30 --concordant --algorithm=blasr --algorithmOptions=--useQuality --maxHits=1 --hitPolicy=random --seed=1 --nproc=8 $line $reference $prefix/$jobid.aln.bam
echo ""
echo "bamtools stats -in $prefix/$jobid.aln.bam"
bamtools stats -in $prefix/$jobid.aln.bam
echo ""
echo "Done"

