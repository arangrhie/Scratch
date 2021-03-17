threads=$1
name=$2
script=$3
args=$4
args_short=`echo $args | awk '{print $1}'`
args_short=`basename $args_short`
script_short=`basename $script`
log=logs/${script_short/.sh/_$args_short.log}

mkdir -p logs
echo "qsub -q phillippy.q -P large_mem -V -pe thread $threads -N $name -cwd -j y -o $log -S /bin/bash $script $args"
qsub -q phillippy.q -P large_mem -V -l mem_free=64G -pe thread $threads -N $name -cwd -j y -o $log -S /bin/bash $script $args

