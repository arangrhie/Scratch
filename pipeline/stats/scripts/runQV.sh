ASM=$1
PREFIX=`basename $ASM`
PREFIX=${PREFIX/.fasta/}
PREFIX=${PREFIX/.fna/}
MAPPING="input/mapping.fofn"
GS=`cat genomesize`

export BWA=/data/projects/phillippy/software/bwa/bwa-0.7.15
export KORENS=/home/korens/devel/utils
export PATH=$KORENS:$BWA:$PATH

ulimit -Su 160000
ulimit -all

PE_LIMIT=5

map=${PREFIX}_map
out_dir=${PREFIX}
mkdir -p $out_dir

#if [ ! -e $PREFIX.busco/success ]; then
#   mkdir -p $PREFIX.busco
#   cd $PREFIX.busco
#   ln -s /data/projects/phillippy/software/BUSCO_v1.22/run.sh
#   ln -s /data/projects/phillippy/software/BUSCO_v1.22/vertebrates
#   ln -s ../$ASM asm.fasta 
#   sh run.sh asm.fasta vertebrates > run.out 2>&1 && touch success
#   cd -
#fi

:<<'END'
if [ ! -e $PREFIX.lumpy.vcf ]; then
   echo "Running lumpy"
   OKTORUN=1
   MISSING=""
   PE=""

   for ill in `cat $MAPPING |grep illumina |awk '{print $1}' |head -n $PE_LIMIT`; do
      illName=`basename $ill`
      BAM="$map/${PREFIX}_${illName}"
      
     # get paired end sizes for each illumina library, output of form mean:344.550742078	stdev:76.3653875471
     readLen=`gunzip -c $ill*_1.fastq.gz 2> /dev/null |head -n 4| tail -n1 | awk '{print length($1)}'`
     if [ x$readLen == "x" ]; then
        readLen=`gunzip -c $ill 2> /dev/null |head -n 4| tail -n1 | awk '{print length($1)}'` 
     fi
     echo "readLen: $readLen"
     pars=`samtools view $BAM.unsorted.bam | tail -n+100000 | /data/projects/phillippy/software/lumpy-sv/scripts/pairend_distro.py -r $readLen -X 4 -N 10000 -o $map/${PREFIX}_${illName}.histo`
     mean=`echo $pars |awk '{print $1}' |sed s/mean://g`
     sd=`echo $pars |awk '{print $NF}' |sed s/stdev://g`

     if [ ! -e $BAM.discordants.bam ]; then
        OKTORUN=0
        MISSING="$MISSING $BAM.discordants.bam"
     fi
     PE="$PE -pe id:sample,bam_file:$BAM.discordants.bam,histo_file:$map/${PREFIX}_${illName}.histo,mean:$mean,stdev:$sd,read_length:$readLen,min_non_overlap:$readLen,discordant_z:5,back_distance:10,weight:1,min_mapping_threshold:20"
   done

#   SR=""
   SR="$SR -sr id:sample,bam_file:$BAM.splitters.bam,back_distance:10,min_mapping_threshold:20,weight:1,min_clip:20"
:<<'END'
   for pac in `cat $MAPPING |grep pacbio |awk '{print $1}'`; do
      pacName=`basename $pac |sed s/.raw.fastq.gz//g |sed s/.fastq.gz//g |sed s/.fastq//g |sed s/.fq.gz//g |sed s/.fq//g`
      BAM="$map/${PREFIX}_${pacName}"
      if [ ! -e $BAM.splitters.bam ]; then
         OKTORUN=0
         MISSING="$MISSING $BAM.splitters.bam"
      else
         SR="$SR -sr id:sample,bam_file:$BAM.splitters.bam,back_distance:10,min_mapping_threshold:20,weight:1,min_clip:20"
      fi
   done
#END
   # run all files
#   if [ $OKTORUN -gt 0 ]; then
#      /data/projects/phillippy/software/lumpy-sv/bin/lumpy -mw 4 -tt 0 $PE > $PREFIX.lumpy.vcf
      /data/projects/phillippy/software/lumpy-sv/bin/lumpy -mw 4 -tt 0 $PE $SR > $PREFIX.lumpy.vcf
#   else
#      echo "Error: $ASM is missing $MISSING inputs running with $SR and $PE"
      #/data/projects/phillippy/software/lumpy-sv/bin/lumpy -mw 4 -tt 0 $PE $SR > $PREFIX.lumpy.vcf
#   fi
fi
END

if [ ! -e ${PREFIX}_FRC.txt ]; then
   echo "Running FRCbam"

   PE=""

   for ill in `cat $MAPPING |grep illumina |awk '{print $1}' |head -n 1`; do
      illName=`basename $ill`
      BAM="$map/${PREFIX}_${illName}"
      echo $BAM
 
      # get insert sizes
      readLen=`gunzip -c $ill*_1.fastq.gz 2> /dev/null |head -n 4| tail -n1 | awk '{print length($1)}'`
      if [ x$readLen == "x" ]; then
         readLen=`gunzip -c $ill 2> /dev/null |head -n 4| tail -n1 | awk '{print length($1)}'`
      fi
      pars=`samtools view $BAM.unsorted.bam | tail -n+100000 | /data/projects/phillippy/software/lumpy-sv/scripts/pairend_distro.py -r $readLen -X 4 -N 10000 -o $map/${PREFIX}_${illName}.histo`
      mean=`echo $pars |awk '{print $1}' |sed s/mean://g`
      sd=`echo $pars |awk '{print $NF}' |sed s/stdev://g`
      MIN=`echo "$mean $sd" |awk '{printf("%d\n", $1-3*$2)}' |awk '{if ($1 < 0) print 0}'`
      MAX=`echo "$mean $sd" |awk '{printf("%d\n", $1+3*$2)}'`

      if [ ! -e $BAM.bam ]; then
         samtools sort -@ 16 -m 2G -o $BAM.bam -T $map/$PREFIX.frc.tmp $BAM.unsorted.bam
      fi

      echo "Using library with insert size $mean $sd ($MIN $MAX)"
      PE="$PE --pe-sam $BAM.bam --pe-max-insert $MAX "
   done
   for ill in `cat $MAPPING |grep illumina |awk '{print $1}' |tail -n 1`; do
      illName=`basename $ill`
      BAM="$map/${PREFIX}_${illName}"
      echo $BAM

      # get insert sizes
      MIN="3000"
      MAX="7000"

      if [ ! -e $BAM.bam ]; then
         samtools sort -@ 16 -m 2G -o $BAM.bam -T $map/$PREFIX.frc.tmp $BAM.unsorted.bam
      fi

      echo "Using mate library with insert size $mean $sd ($MIN $MAX)"
      MP="$MP --mp-sam $BAM.bam --mp-max-insert $MAX "
   done
   echo "/data/projects/phillippy/software/FRC_align/bin/FRC $PE $MP --genome-size $GS --output $out_dir/${PREFIX}"
   /data/projects/phillippy/software/FRC_align/bin/FRC $PE $MP --genome-size $GS --output $out_dir/${PREFIX}
fi
#END

#:<<'END'
# freebayes global params
WEIGHT=0.5
DIPLOID="-p 1 -F $WEIGHT"
if [ -e diploid ]; then
   WEIGHT=0.65
   DIPLOID="-p 2 -F $WEIGHT"
fi

if [ ! -e $PREFIX.bayes.vcf ]; then
   echo "Running vcf"
   if [ ! -e $map/$PREFIX.bam ]; then
      NUM=`cat $MAPPING |grep illumina |wc -l |awk '{print $1}' |head -n $PE_LIMIT`

      MERGE=""

      for ill in `cat $MAPPING |grep illumina |awk '{print $1}' |head -n $PE_LIMIT`; do
            illName=`basename $ill`
            BAM="$map/${PREFIX}_${illName}"
         if [ ! -e $BAM.bam ]; then
            samtools sort -@ 16 -m 2G -o $BAM.bam -T $map/$PREFIX.bayes.tmp $BAM.unsorted.bam
         fi
         if [ ! -e $map/${PREFIX}.header ]; then
            samtools view -H $BAM.unsorted.bam > $map/${PREFIX}.header
         fi

         MERGE="$MERGE $BAM.bam"
      done

      if [ $NUM -gt 1 ]; then
         samtools merge -@16 $map/${PREFIX}.bam $MERGE
         samtools index $map/${PREFIX}.bam
      else
         illName=`cat $MAPPING |grep illumina |head -n 1 |awk '{print $1}'`
         illName=`basename $illName`
         cd $map && ln -s ${PREFIX}_${illName}.bam ${PREFIX}.bam && cd -
	 ln -s ${PREFIX}_${illName}.bam.bai ${PREFIX}.bam.bai
      fi
   fi
   # get SNPs from illumina data
   /data/projects/phillippy/software/freebayes/bin/freebayes -C 2 -0 -O -q 20 -z 0.02 -E 0 -X -u $DIPLOID -b $map/$PREFIX.bam -v $PREFIX.bayes.vcf -f $ASM
fi

if [ ! -e $PREFIX.qv ]; then
   echo "Computing QV (freebayes)"
   # get # snps
   #NUM_SNP=`cat $PREFIX.bayes.vcf |grep -v "#" |awk -v SUM=0 '{if (length($4) == length($5)) { SUM+=length($4); } else if (length($4) < length($5)) { SUM+=length($5)-length($4); } else { SUM+=length($4)-length($5)}} END { print SUM}'`
  NUM_SNP=`cat $PREFIX.bayes.vcf |grep -v "#" | awk -F "\t" '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$8}' | tr ';' ' ' | sed s/AB=//g | awk -v WEIGHT=$WEIGHT '{if ($6 > WEIGHT) print $0}' | awk -v SUM=0 '{if (length($4) == length($5)) { SUM+=length($4); } else if (length($4) < length($5)) { SUM+=length($5)-length($4); } else { SUM+=length($4)-length($5)}} END { print SUM}'`
   # compute QV
   if [ ! -e $PREFIX.numbp ]; then
      # get bases with >3X
      NUM_BP=`samtools depth $map/$PREFIX.bam |awk '{if ($NF >= 3) SUM++; } END { print SUM}'`
      echo "$NUM_BP" > $PREFIX.numbp
   else
      NUM_BP=`cat $PREFIX.numbp`
   fi

   echo "$NUM_SNP" > $PREFIX.numsnp
else
   NUM_BP=`cat $PREFIX.numbp`
   NUM_SNP=`cat $PREFIX.numsnp`
fi
QV=`echo "$NUM_SNP $NUM_BP" | awk '{print (-10*log($1/$2)/log(10))}'`
echo "$QV" > $PREFIX.qv
echo "Num_BP NUM_SNP QV"
echo "$NUM_BP $NUM_SNP $QV"
#END
