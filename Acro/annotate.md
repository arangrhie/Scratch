# Preliminary region annotation of the acrocentric satellites

Annotation will be performed by `annotate.sh` on each assembly, submited via `annotate_job.sh`.

```sh
ls /path/to/*.fa > fa_reforeinted.list

split -l 5 -d -a 2 --numeric-suffixes=10 T2T-chrYv2_fa.list T2T-chrYv2_fa.list.
sh ~/codes/_submit_norm.sh 16 48g annotate annotate_job.sh "T2T-chrYv2_fa.list annotate_v2" "--array=10-30"
```
