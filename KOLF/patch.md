# Patch sequences

## Find corresponding patch sequences
This step is not doable automatically. A newer verkko 2.0 assembly (v0.7) was mapped and examined to the previous v0.6 version,
to find the corresponding distal bits and proximal bits from the rDNA.
`minimap2` was used at first hand to find the best matching scaffolds,
then the assembly graph was further examined along with the raw hi-c interaction maps.

`patch_fa` sequences are under `dima/`, from Dmitry Antipov's verkko 2.0 assembly. If the query sequence contained gaps, to prevent mis-alignments, the sequence were broken at gaps using `break_at_gap.sh`. The following patches were prepared with this script:

* chr14_1
* chr15_1
* chr15_2
* chr22_2

Once the pairing region was identified in a chr_hap pair (`kolf7_patch.list`), for the prior and the newer assembly,
`wfmash --no-split -s100000` was used to get the best matching alignments to v0.6 for each chr_hap of interest using `./patch_to_v0.6.sh`. Each per-patch bam files are named as `chr_hap.n_fix_to_KOLF2.1Jv0.6.bam`. 

At each round, error kmers in the candidate patch sequence was re-evaluated with `./get_error_kmer_pos.sh ${ptch_fa/\.fa/}`.

Merged bam file is `fix_to_KOLF2.1Jv0.6.bam`.

## Adjust target region for patching
Each "fix" patches were evaluated with IGV on each target region.

The error kmers in `*.err.bed` were used to precisely cut to minimize integration of new errors.

Through several rounds of evaluation, `patch_candidate.txt` was generated, which contains the following:
```
reference	query	strand	query_chr_hap	[query	strand	query_chr_hap]
```
If a second [query	strand	query_chr_hap] was provided, it is assumed to place a 1Mbp gap between the first and the second querey sequences.

Patch sequences were then generated in a VCF format, to contain the original sequence in REF field, and the replacement sequence in ALT field. Final VCF `KOLF2.1Jv0.6.patch.vcf.gz` was then used to generated a new candidate consensus `patched.fa`.

```
bcftools consensus -c $sp_ver.to_patched.chain -f $sp_ver.fa -HA $sp_ver.patch.vcf.gz > patched.fa
```

Patched target regions +- 50 kbp were extracted and lifted over to the new patched.fa coordinates to extract the final patches +- 50 kbps and merged. Telomeric patches were manually edited in `patch_50k.lifted.mrg.bed` to accomodate unlifted regions.

Each region was then used to extract patched sequences and re-aligned to v0.6 assembly to confirm patched sequences using `wfmash --no-split -ad -t$cpus ${sp_ver}.$chr_hap.fa \
    patch.$i.fa -s20000 -p95 --one-to-one > $out.sam`.

Final bam file containing all patched sequences aligned to v0.6 is named `patch_to_KOLF2.1Jv0.6.bam`.

These steps above are performed in the following script:
```shell
./make_patch.sh
```

An IGV session for the patch alignment is available as `igv_patch_to_KOLF2.1Jv0.6.xml` using paths under CARD_AUX.

## Merge with SNV edits and apply correction
Under `/data/CARD_AUX/KOLF/polish/consensus`
* SV: ../patch/KOLF2.1Jv0.6.patch.vcf.gz
* SNV: ../variant/snv_candidates.merfin-loose.vcf.gz
```shell
./apply_correction.sh
```

Changes to v1.0 candidate:
* Revert chrX SV correction - more error kmers than SNV based correction.
* Fixing `chr13_2` - qry region was invalid; fixing `haplotype1-0000015:0-355887` to `haplotype1-0000015:1-355887`
* Intermediate assembly v0.9.1 applying these changes before haplotype switching