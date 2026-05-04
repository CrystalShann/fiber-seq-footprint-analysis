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
  - Reads combined FIRE nucleosome intervals from `/project/spott/cshan/fiber-seq/FIRE_nuc_by_chr_combined_sample/parsed_bed_files/sorted_indexed`
  - Builds clipped per-read BED12 tracks for m6A plus footprint size bins `10-30`, `40-60`, `60-80`, and `140-160`
  - Uses Hi-C from `/project/spott/cshan/annotations/hic/ENCFF216QQM.mcool` for both small and large regions
  - Converts `/project/spott/cshan/annotations/hic/ENCFF156ECM.bedpe.gz` to a merged BED boundary file under each region output folder and uses that BED for small regions
  - Uses CAGE peaks from `/project/spott/cshan/fiber-seq/results/TADs_m6a_footprint_overlay/inputs/cage_lcl_extended.bed.gz`
  - Uses FIRE combined nucleosome coverage bigWigs from `/project/spott/cshan/fiber-seq/FIRE_nuc_by_chr_combined_sample/bigwig_files`
  - Discovers PolII footprint coverage tracks from `/project/spott/cshan/fiber-seq/results/PolII/FIRE_combined_footprints/joint_trained_tracks`
    by chromosome and footprint range parsed from filenames like `combined_chr1_100-120bp_fp_cov.bw`
  - Writes a region-specific `pyGenomeTracks` `.ini`
- `run_region_overlay.sh`
  - Runs both steps
  - Renders the final PNG and PDF with `pyGenomeTracks`

Example:

```bash
bash /project/spott/cshan/fiber-seq/code/TADs_m6a_footprint_overlay_code/run_region_overlay.sh \
  --sample AL10_bc2178_19130 \
  --region chr1:2320000-2340000

sbatch /project/spott/cshan/fiber-seq/code/TADs_m6a_footprint_overlay_code/run_region_overlay.sh \
  --sample AL60_bc2152_18519 \
  --region chr1:1724000-1726000
```

Outputs land in:

```text
/project/spott/cshan/fiber-seq/results/TADs_m6a_footprint_overlay/<sample>/<chrom>/<chrom>_<start>_<end>/
```

Track order (top to bottom):

1. x-axis
2. Hi-C triangle (ENCFF216QQM mcool; always shown)
3. TAD boundaries (converted BED from ENCFF156ECM.bedpe.gz for small regions) or TAD domains (Arrowhead, GM12878)
4. CTCF bigWig
5. DNase bigWig
6. Nucleosome coverage bigWig (FIRE combined)
7. FANTOM5 CAGE
8. CpG methylation — `fire_CpG/<sample>/<sample>_CPG.combined.bed.gz` rendered as a bedgraph (scale 0–100); included automatically when the file exists
9. PolII footprint coverage bigWigs discovered from `results/PolII/FIRE_combined_footprints/joint_trained_tracks` for the requested chromosome; footprint ranges are parsed from filenames like `combined_chr1_100-120bp_fp_cov.bw`
10. m6A fraction — aggregate from `footprintR` on the BAM (purple)
11. 5mC fraction — aggregate from `footprintR` on the BAM (blue); included when the BAM carries 5mC tags
12. Per-read overlay — baseline gray reads with m6A (purple) and footprint size bins overlaid

Notes:

- The CpG track uses the pb-CpG-tools combined BED file directly (column 4 = model score, 0–100) instead of the derived bigWig.
- The 5mC aggregate track is computed from the BAM using the same `footprintR::readModBam` approach as m6A, with `modbase_code = "m"`.
- For regions smaller than 500 kb, the track builder uses converted BED TAD boundaries while keeping Hi-C enabled.
