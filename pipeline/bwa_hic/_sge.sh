threads=$1
name=$2
script=$3
args=$4
args_simple=`basename $args`
log=logs/${3/.sh/_$args_simple.log}

mkdir -p logs

echo "qsub -q phillippy.q -pe thread $threads -N $name -cwd -o $log -e $log -S /bin/bash $script $args"
qsub -q phillippy.q -pe thread $threads -N $name -cwd -o $log -e $log -S /bin/bash $script $args

