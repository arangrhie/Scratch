# Find missing telomere

Find and retrieve sequences with telomeres intact
```
cut -f1 telomere.bed | sort -u > telomere.list
samtools faidx -@12 -r telomere.list assembly.fasta > kolf7_dima.telo.fa
```

Map on v1.0
```shell
module load wfmash/0.14.0
module load samtools

out=telo_to_v1.0
wfmash --no-split -ad -t24 KOLF2.1Jv1.0.fa \
    kolf7_dima.telo.fa -s100000 -p95 --one-to-one > $out.sam
samtools sort -@24 -T $out.tmp -O bam $out.sam > $out.bam
samtools index $out.bam

# For the 1 node < 100 kb
awk '$2<100000 {print $1}' kolf7_dima.telo.fa.fai > kolf7_dima.telo_under100k.list
samtools faidx -r kolf7_dima.telo_under100k.list kolf7_dima.telo.fa > kolf7_dima.telo_under100k.fa

out=telo_short_to_v1.0
wfmash --no-split -ad -t24 KOLF2.1Jv1.0.fa \
    kolf7_dima.telo_under100k.fa -s50000 -p95 --one-to-one > $out.sam
samtools sort -@24 -T $out.tmp -O bam $out.sam > $out.bam
samtools index $out.bam
```

## Collect qry cov, ref cov, and if it's reaching the p-tel or q-tel

```sh
## Collect rough p|q tel presence
cut -f1-12 telo_to_v1.0.paf |\
  awk -v OFS="\t" '{ ptel=""; qtel=""; \
       { if ($8 < 10000) { ptel="p" }; } \
       { if ( ($7 - $9) < 10000) { qtel = "q" }; } \
       { print $0, 100*$(NF-2)/$2, 100*$(NF-1)/$7, ptel "|" qtel } }'

## Collect those that are covered nearly t2t
cut -f1-12 telo_to_v1.0.paf |\
  awk -v OFS="\t" '{ ptel=""; qtel=""; \
       { if ($8 < 10000) { ptel="p" }; } \
       { if ( ($7 - $9) < 10000) { qtel = "q" }; } \
       { print $0, 100*$(NF-2)/$2, 100*$(NF-1)/$7, ptel "|" qtel } }' |\
  awk '{if ($(NF-2) > 90 && $(NF-1) > 90 ) print $0}'

## Collect the rest
cut -f1-12 telo_to_v1.0.paf |\
  awk -v OFS="\t" '{ ptel=""; qtel=""; \
       { if ($8 < 10000) { ptel="p" }; } \
       { if ( ($7 - $9) < 10000) { qtel = "q" }; } \
       { print $0, 100*$(NF-2)/$2, 100*$(NF-1)/$7, ptel "|" qtel } }' |\
  awk -v OFS="\t" '{ if ($(NF-2) <= 90 || $(NF-1) <= 90 ) print $0}'
```

`unassigned-0000410` is the only one not assigned for the telomeres. Makes it a potential target for the missing telomere end!!

* map unassigned-0000410 to chr21_2
![alt text](image.png)



```sh
module load samtools/1.21
module load minimap2/2.28

samtools faidx KOLF2.1Jv1.0.fa chr21_2 > chr21_2.fa
samtools faidx KOLF2.1Jv1.0.fa chr22_2:1-400000 > chr22_2.fa
samtools faidx -i assembly.fasta unassigned-0000410 > telo_patch_21_2p.fa

minimap2 -xasm5 --MD -c -t 24 chr21_2.fa telo_patch_21_2p.fa > telo_patch_21_2p.to.chr21_2.paf
```
![alt text](image-1.png)
* Confirmed it's poorly mapping to chr21_2.
* Check how identical the two contigs are 
```sh
minimap2 -xasm10 --MD -c -t 24 chr22_2.fa telo_patch_21_2p.fa > telo_patch_21_2p.to.chr22_2.paf
cut -f1-12 telo_patch_21_2p.to.chr22_2.paf | sort -k3,3n
unassigned-0000410/rc   393475  164     11390   +       chr22_2:1-400000        400000  31      11273   11191   11269   24
unassigned-0000410/rc   393475  930     11390   +       chr22_2:1-400000        400000  31      10505   10428   10496   0
unassigned-0000410/rc   393475  1452    11744   +       chr22_2:1-400000        400000  2170    12342   10123   10311   0
unassigned-0000410/rc   393475  1453    11503   +       chr22_2:1-400000        400000  2320    12371   10000   10068   0
unassigned-0000410/rc   393475  14621   354640  +       chr22_2:1-400000        400000  12769   209451  193381  342229  60
unassigned-0000410/rc   393475  354844  377198  +       chr22_2:1-400000        400000  281530  317903  22005   36381   60

sed 's/chr22_2:1-400000/chr22_2/g' telo_patch_21_2p.to.chr22_2.paf > telo_patch_21_2p.to.chr22_2.adjust.paf

## Also check chr14_2
samtools faidx KOLF2.1Jv1.0.fa chr14_2:1-400000 > chr14_2.fa
minimap2 -xasm10 --MD -c -t 24 chr14_2.fa telo_patch_21_2p.fa > telo_patch_21_2p.to.chr14_2.paf
cut -f1-12 telo_patch_21_2p.to.chr14_2.paf | sort -k7,7n

unassigned-0000410/rc   393475  1642    303885  +       chr14_2:1-400000        400000  1748    309178  299823  308704  60
unassigned-0000410/rc   393475  81      1429    +       chr14_2:1-400000        400000  10280   11672   1306    1414    60
unassigned-0000410/rc   393475  81      883     +       chr14_2:1-400000        400000  10802   11672   773     877     0

sed 's/chr14_2:1-400000/chr14_2/g' telo_patch_21_2p.to.chr14_2.paf > telo_patch_21_2p.to.chr14_2.adjust.paf

minimap2 -xasm10 --MD -a -t 24 chr14_2.fa telo_patch_21_2p.fa > telo_patch_21_2p.to.chr14_2.sam
sed 's/chr14_2:1-400000/chr14_2/g' telo_patch_21_2p.to.chr14_2.sam > telo_patch_21_2p.to.chr14_2.adjust.sam
```
![alt text](image-2.png)
Maps nicely on chr14_2, possibly this actually belongs here.
However, given the diverged sequences towards the proximal site of the satellite, sequence divergence happens. Could be just two very similar distal sequences.

## Proceeding in `../patch/`
Generate `chr21_2.new.fa` that will be replaced with the chr21_2 in v1.0, and bumped up to v1.1

```sh
cd ../patch

ln -s ../telomere/telo_patch_21_2p.fa

./get_error_kmer_pos.sh telo_patch_21_2p
cat telo_patch_21_2p.err.bed 
# unassigned-0000410/rc   15      109

# 10kb = 100 * 100
for i in $(seq 1 100);
do
  if [[ $(($i % 100)) -eq 0 ]]; then
    echo $i;
  fi
  cat gap.100N.seq  >> gap.seq;
done
echo "N" >> gap.seq
mv gap.seq gap.10001.seq

cd ../
mkdir patch_v1.0_chr21_2
cd patch_v1.0_chr21_2

ln -s ../telomere/telo_patch_21_2p.fa
ln -s ../telomere/telo_patch_21_2p.fa.fai
ln -s ../patch/gap.10001.seq

ln -s ../telomere/chr21_2.fa
ln -s ../telomere/chr21_2.fa.fai

tmp_seq=tmp.seq

samtools faidx telo_patch_21_2p.fa unassigned-0000410/rc:110-393475 | awk 'NR>1' - > $tmp_seq # 393366 bp

# Add gap
cat gap.10001.seq >> $tmp_seq

# Add the rest of the chr21_2
awk 'NR>1' chr21_2.fa >> $tmp_seq

echo ">chr21_2" > chr21_2.tmp.fa
cat $tmp_seq | tr -d '\n' >> chr21_2.tmp.fa

samtools faidx chr21_2.tmp.fa
samtools faidx -n60 chr21_2.tmp.fa chr21_2 > chr21_2.new.fa

ln -s /data/CARD_AUX/KOLF/polish/assemblies/v1.0/v1.0.hap2.fa
ln -s /data/CARD_AUX/KOLF/polish/assemblies/v1.0/v1.0.hap2.fa.fai

# collect everything from chr1 to chr20
for i in $(seq 1 20)
do
  echo "chr${i}_2" >> seq_tmp.1
done

# and for chr22 and Y
echo "chr22_2" > seq_tmp.2
echo "chrY"   >> seq_tmp.2

samtools faidx -r seq_tmp.1 -@12 v1.0.hap2.fa >  v1.1.hap2.fa
cat chr21_2.new.fa >> v1.1.hap2.fa
samtools faidx -r seq_tmp.2 -@12 v1.0.hap2.fa >> v1.1.hap2.fa
samtools faidx -@12 v1.1.hap2.fa

# confirm the length
cat v1.1.hap2.fa.fai
## chr21_2 44395731        2828148257      60      61
## v1.0
## chr21_2 43992364        2828148257      60      61

```