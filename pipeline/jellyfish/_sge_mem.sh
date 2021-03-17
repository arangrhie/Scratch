threads=$1
mem=$2
name=$3
script=$4
args=$5
script_basename=`basename $script`
args_simple=`basename $args`
log=logs/${script_basename/.sh/_$args_simple.log}

mkdir -p logs

echo "qsub -V -q phillippy.q -P large_mem -pe thread $threads -l mem_free=$mem -N $name -cwd -o $log -e $log -S /bin/bash $script $args"
qsub -V -q phillippy.q -P large_mem -pe thread $threads -l mem_free=$mem -N $name -cwd -o $log -e $log -S /bin/bash $script $args

