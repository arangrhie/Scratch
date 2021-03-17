threads=$1
name=$2
script=$3
args=$4
script_basename=`basename $script`
args_simple=`basename $args`
log=logs/${script_basename/.sh/_$args_simple.log}

mkdir -p logs

echo "qsub -V -q phillippy.q -P large_mem -pe thread $threads -l mem_free=45g -N $name -cwd -o $log -e $log -S /bin/bash $script $args"
qsub -V -q phillippy.q -P large_mem -pe thread $threads -l mem_free=45g -N $name -cwd -o $log -e $log -S /bin/bash $script $args

