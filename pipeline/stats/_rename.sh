rename=$1
to=$2

for file in $(ls -d ./${rename}*); do
	echo "mv $file ${file/$rename/$to}"
	mv $file ${file/$rename/$to}
done
