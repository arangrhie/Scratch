threads=$1

#/data/projects/phillippy/software/samtools-1.3/bin/samtools sort -n -o pac_fcn.10x.sam -m2g -T pac_fcn.10x.tmp -@ $threads -O sam pos_sorted_bam.bam

#java -jar -Xmx1g /home/rhiea/codes/tenXGtoArcs.jar pac_fcn.10x.sam | samtools view -hb -@$threads - > pac_fcn.10x.arcs.bam

/data/projects/phillippy/software/arcs/bin/arcs -f pac_fcn_p.fasta -a bam_list.txt -b pac_fcn_p_10x_arcs

/data/projects/phillippy/software/arcs/makeTSVfile.py pac_fcn_p_10x_arcs_original.gv pac_fcn_p_10x_arcs.tigpair_checkpoint.tsv pac_fcn_p.fasta

touch empty.fof
/data/projects/phillippy/software/links_v1.8.5/LINKS -f pac_fcn_p.fasta -s empty.fof -b pac_fcn_p_10x_arcs -k 20
