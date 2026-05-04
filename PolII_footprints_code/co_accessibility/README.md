# Co-accessibility around high and low pausing genes

This folder contains an R Markdown workflow that adapts Kevin's co-accessibility analysis to high and low pausing gene sets.

## Files
- coaccess_fire_cres_high_low.Rmd: main analysis notebook
- build_coaccess_inputs.R: builds `genomic.annots.hg38.rds`, high/low pausing gene lists, and a ranked PI table from the CAGE-TSS pausing-index table
- extract_region_fire_result.sh: local region extraction helper used by the Rmd

## Required inputs
- A gene annotation RDS with `genomic.annots$genes` as a GRanges and a `gene_name` column.
- High and low pausing gene lists (one gene symbol per line).
- ENCODE cCREs BED (SCREEN v4).
- Fiber-seq FIRE outputs for the sample (peaks and fire-elements).
- BAM file for the sample.

## Build co-accessibility inputs
Run this once before rendering the Rmd:

```bash
Rscript /project/spott/cshan/fiber-seq/code/PolII_footprints_code/co_accessibility/build_coaccess_inputs.R 500
```

The optional numeric argument is the number of genes to write in each list. The script uses:

- Pausing index table: `/project/spott/cshan/fiber-seq/results/PolII/annotations/pausing_index_principal_with_CAGE_TSS.csv`
- ENCODE cCRE BED: `/project/spott/cshan/annotations/GRCh38-cCREs.bed`
- Output directory: `/project/spott/cshan/fiber-seq/results/co_accessibility/inputs`

The gene annotation RDS stores `genomic.annots$genes` as a GRanges with `gene_name`, `gene_id`, `PI`, and a `tss` metadata column that prefers `cage_tss` and falls back to the table's original `tss`.

## Parameters you must set in the Rmd
- `sample_name`, `bam_file`, `fire_dir`, `sample_fire_dir`, `out_dir`
- `genomic_annots_rds`, `high_gene_list_file`, `low_gene_list_file`, `cre_bed`
- `selected_chr`, `tss_flank`, `min_dist`, `max_dist`, `pseudocount`, `ncore`

## Outputs
For each gene set (high or low) and each chromosome:
- `*_CREs_around_genes_coaccess_res.rds`: full list of per-gene co-access results
- `*_CREs_around_genes_coaccess_stat_df.rds`: flat per-pair stats with counts, OR, and p-value

## Notes
- The workflow now includes local copies of `extract_region_fire_result`, `extract_fire_read_info`, and `test_coaccess_regions`, so it does not source Kevin's `process_fiberseq_data.R`.
- Run on a single chromosome first to validate paths and runtime before scaling to all chromosomes.
