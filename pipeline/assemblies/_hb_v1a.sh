sort -k2 -nr len/pac_nano_canu_arrow_arrow_pilon_10x_spnv_hap1_qm2_arrow_bio_tgh_all_salsa1_pbjl_bio_tgh_all_rename.fasta.len | head -n 500 | awk '{print $1}' > map/hb_v1a.fasta.list

