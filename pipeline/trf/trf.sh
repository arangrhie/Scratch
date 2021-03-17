### trf.sh

ASM=$1
PREFIX=${ASM/.fasta/}
PREFIX=${PREFIX/.fna/}
PREFIX=${PREFIX/.fa/}

LEN="$PREFIX.lens"

SCRIPT_PATH=/data/projects/phillippy/software/PilonGrid
JAVA_PATH=$SCRIPT_PATH:.

if [ ! -e $LEN ]; then
   java -cp $SCRIPT_PATH:. SizeFasta $ASM > $LEN
fi
NUM_CTG=`wc -l $LEN |awk '{print $1}'`
echo "NUM_CTG: $NUM_CTG"

mkdir -p logs
qsub -V -pe thread 1 -q phillippy.q -tc 400  -l mem_free=2G -t 1-$NUM_CTG -cwd -N "trf.${PREFIX}" -j y -o logs/$PREFIX.\$TASK_ID.log -S /bin/bash trf_parallel.sh $PREFIX
