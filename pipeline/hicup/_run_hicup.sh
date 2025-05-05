./scripts/hicup_digester --re1 A^GATCT,BglII --genome pac_falcon_pa_pilon_ponly_hicup_GATC --outdir pac_falcon_pa_pilon_ponly_hicup_GATC pac_falcon_pa_pilon_ponly/pacbio_falcon_pa_pilon_ponly.fasta
./scripts/hicup_digester --re1 A^GATCT,BglII:A^GAATCT,arima1:A^GAGTCT,arima2:A^GATTCT,arima3:A^GACTCT,arima4 --genome pac_falcon_pa_pilon_ponly_hicup_GATC_GANTC --outdir pac_falcon_pa_pilon_ponly_hicup_GATC_GANTC pac_falcon_pa_pilon_ponly/pacbio_falcon_pa_pilon_ponly.fasta
./_sge.sh 24 hicup_dovetail8 _hicup.sh dovetail8
./_sge.sh 24 hicup_dovetail9 _hicup.sh dovetail9
./_sge.sh 24 hicup_phase _hicup.sh phase
./_sge.sh 24 hicup_arima _hicup.sh arima
