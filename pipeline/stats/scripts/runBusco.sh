#for ASM in $(cat input/assembly.list); do
ASM=$1
cpus="-c $2"
PREFIX=`basename $ASM`
PREFIX=${PREFIX/.fasta/}

ulimit -Su 160000
ulimit -all

echo "Running BUSCO on $ASM in $PREFIX.busco"
if [ ! -e $PREFIX.busco/success ]; then
  mkdir -p $PREFIX.busco
  cd $PREFIX.busco
  if [ ! -e run.sh ] ; then
    ln -s /data/projects/phillippy/software/BUSCO_v1.22/run.sh
    ln -s /data/projects/phillippy/software/BUSCO_v1.22/vertebrates
  fi
  echo "sh run.sh $ASM vertebrates $cpus > run.out 2>&1 && touch success"
  sh run.sh $ASM vertebrates "$cpus" > run.out 2>&1 && touch success
  cd -
fi
#done
