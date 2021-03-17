assembly=$1
assembly=${assembly/.fasta/}
assembly=`basename $assembly`

resite="GANTC"	# arima GATC and GANTC

cd restriction_sites
if [ ! -e "$assembly.fasta" ]; then
	ln -s ../references/$assembly.fasta
fi
if [ ! -e "$assembly.$resite.txt" ] ; then
	echo "Generating $assembly.$resite.bed"
	java -jar -Xmx4g /home/rhiea/codes/fastaFindSequence.jar $resite $assembly.fasta > $assembly.$resite.bed
	java -jar -Xmx1g /home/rhiea/codes/juicerREbedToTxt.jar $assembly.$resite.bed > $assembly.$resite.txt
fi
if [ ! -e "$assembly.chr.sizes" ] ; then
	echo "Generating $assembly.chr.sizes"
	java -jar -Xmx1g /home/rhiea/codes/fastaContigSize.jar $assembly.fasta
	awk '{print $1" "$2}' $assembly.fasta.len > $assembly.chr.sizes
fi
cd ../

mkdir -p $assembly
cd $assembly
if ! [ -e references ]; then
	ln -s ../references
	ln -s ../restriction_sites
	ln -s ../fastq
fi
cd ..

if ! [ -e $assembly/aligned/merged_nodups.txt ]; then
#echo "/data/projects/phillippy/software/juicer/scripts/juicer.sh -g hbs2 -z $PWD/references/$assembly.fasta -y $PWD/restriction_sites/$assembly.$resite.txt -q phillippy.q -l phillippy.q -D /data/projects/phillippy/software/juicer/ -d $PWD/$assembly -p $PWD/restriction_sites/$assembly.chr.sizes"
#/data/projects/phillippy/software/juicer/scripts/juicer.sh -g hbs2 -z $PWD/references/$assembly.fasta -y $PWD/restriction_sites/$assembly.$resite.txt -q phillippy.q -l phillippy.q -D /data/projects/phillippy/software/juicer/ -d $PWD/$assembly -p $PWD/restriction_sites/$assembly.chr.sizes
echo "/data/projects/phillippy/software/juicer/scripts/juicer.sh -S early -g hbs2 -z $PWD/references/$assembly.fasta -y $PWD/restriction_sites/$assembly.$resite.txt -q phillippy.q -l phillippy.q -D /data/projects/phillippy/software/juicer/ -d $PWD/$assembly -p $PWD/restriction_sites/$assembly.chr.sizes"
/data/projects/phillippy/software/juicer/scripts/juicer.sh -S early -g hbs2 -z $PWD/references/$assembly.fasta -y $PWD/restriction_sites/$assembly.$resite.txt -q phillippy.q -l phillippy.q -D /data/projects/phillippy/software/juicer/ -d $PWD/$assembly -p $PWD/restriction_sites/$assembly.chr.sizes
fi

if [ -e $assembly/aligned/merged_nodups.txt ]; then
	cd $assembly
	mkdir -p scaff
	cd scaff
	if ! [ -e asm.fasta ]; then
		echo "Create symlinks"
		ln -s ../aligned/merged_nodups.txt
		ln -s ../references/$assembly.fasta asm.fasta
		ln -s ../../_3ddna.sh
	fi
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
	   echo "Error, python didn't work!"
	   exit
	fi

	sh /home/rhiea/codes/_sge.sh 64 juicer3ddna _3ddna.sh
fi
