threads=$1
name=$2
script=$3
args=$4
args_short=`basename $4`
log=logs/${3/.sh/_$args_short.log}

mkdir -p logs

echo "qsub -q phillippy.q -P large_mem -pe thread $threads -N $name -cwd -o $log -e $log -S /bin/bash $script $args"
qsub -q phillippy.q -P large_mem -pe thread $threads -N $name -cwd -o $log -e $log -S /bin/bash $script $args

