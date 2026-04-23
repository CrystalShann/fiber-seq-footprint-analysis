# Region m6A + Footprint Overlay

This folder adds a sample-aware region workflow that keeps the existing `pyGenomeTracks` stack and appends a per-read overlay panel.

Files:

- `extract_region_m6a_summary.R`
  - Looks up the sample in `LCL_sample_metatable_merged_samples_31samples.csv`
  - Reads the sample BAM
  - Computes region-level m6A using `footprintR::readModBam(..., modbase_code = "a", level = "summary")`
  - Writes `m6A_position_calls.tsv.gz`, `m6A_fraction.bedgraph`, optional `m6A_fraction.bw`, and a small PDF summary
  - Also computes region-level 5mC using `footprintR::readModBam(..., modbase_code = "m", level = "summary")`
  - Writes `5mC_fraction.bedgraph` and optional `5mC_fraction.bw`
- `build_region_overlay_tracks.py`
  - Reads the sample-specific FIRE `m6a_by_chr` BED12
  - Reads the sample-specific FiberHMM `joint_trained_results/*_chr*_fp.bed.gz`
  - Builds clipped per-read BED12 tracks for m6A plus footprint size bins `10-30`, `40-60`, `60-80`, and `140-160`
  - Writes a region-specific `pyGenomeTracks` `.ini`
- `run_region_overlay.sh`
  - Runs both steps
  - Renders the final PNG and PDF with `pyGenomeTracks`

Example:

```bash
bash /project/spott/cshan/fiber-seq/code/region_m6a_footprint_overlay/run_region_overlay.sh \
  --sample AL10_bc2178_19130 \
  --region chr1:2320000-2340000
```

Outputs land in:

```text
/project/spott/cshan/fiber-seq/results/region_m6a_footprint_overlay/<sample>/<chrom>/<chrom>_<start>_<end>/
```

Track order (top to bottom):

1. x-axis
2. Hi-C triangle (GM12878, 5 kb)
3. TADs (Arrowhead, GM12878)
4. CTCF bigWig
5. DNase bigWig
6. FANTOM5 CAGE
7. CpG methylation — `fire_CpG/<sample>/<sample>_CPG.combined.bed.gz` rendered as a bedgraph (scale 0–100); included automatically when the file exists
8. Footprint coverage bigWigs (10-30, 40-60, 60-80, 140-160 bp)
9. m6A fraction — aggregate from `footprintR` on the BAM (purple)
10. 5mC fraction — aggregate from `footprintR` on the BAM (blue); included when the BAM carries 5mC tags
11. Per-read overlay — baseline gray reads with m6A (purple) and footprint size bins overlaid

Notes:

- The CpG track uses the pb-CpG-tools combined BED file directly (column 4 = model score, 0–100) instead of the derived bigWig.
- The 5mC aggregate track is computed from the BAM using the same `footprintR::readModBam` approach as m6A, with `modbase_code = "m"`.
