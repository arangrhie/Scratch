export BWA=/data/projects/phillippy/software/bwa/bwa-0.7.15
export PATH=$BWA:$PATH

threads=$1
samtools view -h dovetail_0799_to_pac_falcon.r1.bam | awk '{if ($1 !~ /^@/) {$2 += 65; print $0} else {print $0}}' - > dovetail_0799_to_pac_falcon.fix.sam
samtools view dovetail_0799_to_pac_falcon.r2.bam | awk '{if ($1 !~ /^@/) {$2 += 129; print $0} else {print $0}}' - >> dovetail_0799_to_pac_falcon.fix.sam
samtools sort -m2g -@$threads -n -O bam -o dovetail_0799_to_pac_falcon.fix.bam -T dovetail_0799_to_pac_falcon.fix.tmp dovetail_0799_to_pac_falcon.fix.sam
samtools fixmate dovetail_0799_to_pac_falcon.fix.bam dovetail_0799_to_pac_falcon.fix.fixmate.bam
samtools sort -m2g -@$threads -O bam -o dovetail_0799_to_pac_falcon.fix.fixmate.sort.bam -T dovetail_0799_to_pac_falcon.fix.tmp dovetail_0799_to_pac_falcon.fix.fixmate.bam
samtools stats dovetail_0799_to_pac_falcon.fix.fixmate.sort.bam > dovetail_0799_to_pac_falcon.fix.fixmate.sort.stats
