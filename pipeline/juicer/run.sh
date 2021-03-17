export TMPDIR=`pwd`/tmp
mkdir -p $TMPDIR

module load Bio/lastz/1.04.00
module load util/coreutils/8.27
export PATH=/data/projects/phillippy/software/gawk-4.0.2/bin:$PATH
# test module is ok and python version is ok
module load libs/scipy/0.18.1 > check 2>&1

NUM=`cat check |wc -l |awk '{print $1}'`
if [ $NUM -ne 0 ]; then
   echo "Error cannot run, python failing to load scipy!"
   exit
fi

export PYTHONPATH=/data/projects/phillippy/software/scipy-0.18.1/lib/python2.7/site-packages/:$PYTHONPATH
which python
python /data/projects/phillippy/software/3d-dna/checkModule.py
if [ $? -ne 0 ]; then
   echo "error python didn't work!"
   exit
fi

bash /data/projects/phillippy/software/3d-dna/run-asm-pipeline.sh -m haploid -i 15000 -r 2  asm.fasta merged_nodups.txt
