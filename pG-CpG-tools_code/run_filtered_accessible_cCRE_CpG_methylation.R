required_pkgs <- c(
  "data.table",
  "dplyr",
  "tidyr",
  "tibble",
  "stringr",
  "GenomicRanges",
  "IRanges",
  "GenomeInfoDb",
  "SummarizedExperiment",
  "BiocParallel",
  "Rsamtools",
  "footprintR",
  "ggplot2",
  "patchwork",
  "scales",
  "signal"
)

missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_pkgs) > 0) {
  stop(
    sprintf(
      "Missing required packages: %s",
      paste(missing_pkgs, collapse = ", ")
    ),
    call. = FALSE
  )
}

invisible(lapply(required_pkgs, library, character.only = TRUE))

# ── Variables ─────────────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
default_sample <- ifelse(length(args) >= 1, args[[1]], "AL10_bc2178_19130")
default_chrom  <- ifelse(length(args) >= 2, args[[2]], "chr1")

ccre_output_dir <- "/project/spott/cshan/fiber-seq/results/cCRE_summary/filtered_DNase_cCRE_mCpG"
dir.create(ccre_output_dir, recursive = TRUE, showWarnings = FALSE)

pbcpg_bed_path <- file.path(
  "/project/spott/cshan/fiber-seq/results/fire_CpG",
  default_sample,
  sprintf("%s_CPG.combined.bed.gz", default_sample)
)
filtered_ccre_path <- file.path(ccre_output_dir, "DNase_overlapping_cCREs.tsv.gz")
filtered_ccre_bed_path <- file.path(ccre_output_dir, "DNase_overlapping_cCREs.bed")
class_summary_tsv_path <- file.path(ccre_output_dir, "DNase_overlapping_cCREs_class_summary.tsv")
class_summary_gz_path  <- file.path(ccre_output_dir, "DNase_overlapping_cCREs_class_summary.tsv.gz")
dnase_path <- "/project/spott/cshan/annotations/DNase_ENCFF073ORT.bed"
ccre_path  <- "/project/spott/cshan/annotations/GRCh38-cCREs.bed"

# ── Helpers ───────────────────────────────────────────────────────────────────
assert_file_exists <- function(path, label = "file") {
  if (!file.exists(path)) stop(sprintf("Missing %s: %s", label, path), call. = FALSE)
  invisible(path)
}

write_gz_tsv <- function(x, path) {
  data.table::fwrite(x, file = path, sep = "\t")
  invisible(path)
}

summarize_region_pbcpg_5mc <- function(region_row, pbcpg_bed) {
  hits <- pbcpg_bed[
    chrom == region_row$chromosome &
      begin >= region_row$start0 &
      end <= region_row$end0
  ]
  modified_calls   <- as.integer(round(sum(hits$est_mod_count)))
  unmodified_calls <- as.integer(round(sum(hits$est_unmod_count)))
  total_calls      <- modified_calls + unmodified_calls
  data.table::data.table(
    region_class = region_row$region_class,
    chromosome = region_row$chromosome,
    start = region_row$start0 + 1L,
    end = region_row$end0,
    region_width_bp = region_row$end0 - region_row$start0,
    n_cpg_sites = nrow(hits),
    modified_calls = modified_calls,
    unmodified_calls = unmodified_calls,
    total_calls = total_calls,
    fraction_modified = data.table::fifelse(
      total_calls > 0L, modified_calls / total_calls, NA_real_
    )
  )
}

# ── Load data ─────────────────────────────────────────────────────────────────
assert_file_exists(pbcpg_bed_path, "pb-CpG-tools combined BED")
if (!file.exists(filtered_ccre_path)) {
  assert_file_exists(dnase_path, "DNase BED")
  assert_file_exists(ccre_path,  "cCRE BED")
  
  ccre_dt <- data.table::fread(
    ccre_path,
    header = FALSE,
    col.names = c("chromosome", "start0", "end0", "ccre_id", "ccre_element_id", "ccre_class")
  )
  dnase_dt <- data.table::fread(
    dnase_path,
    header = FALSE,
    select = 1:3,
    col.names = c("chromosome", "start0", "end0")
  )
  
  ccre_gr <- GenomicRanges::GRanges(
    seqnames = ccre_dt$chromosome,
    ranges   = IRanges::IRanges(start = ccre_dt$start0 + 1L, end = ccre_dt$end0)
  )
  dnase_gr <- GenomicRanges::GRanges(
    seqnames = dnase_dt$chromosome,
    ranges   = IRanges::IRanges(start = dnase_dt$start0 + 1L, end = dnase_dt$end0)
  )
  
  hits <- GenomicRanges::findOverlaps(ccre_gr, dnase_gr, ignore.strand = TRUE)
  ccre_idx <- sort(unique(S4Vectors::queryHits(hits)))
  filtered_ccre_dt_build <- ccre_dt[ccre_idx]
  
  summary_dt <- filtered_ccre_dt_build[, .N, by = ccre_class][order(-N)]
  summary_dt[, fraction := N / sum(N)]
  
  data.table::fwrite(
    filtered_ccre_dt_build,
    filtered_ccre_bed_path,
    sep = "\t",
    col.names = FALSE
  )
  data.table::fwrite(
    summary_dt,
    class_summary_tsv_path,
    sep = "\t"
  )
  data.table::fwrite(
    summary_dt,
    class_summary_gz_path,
    sep = "\t"
  )
  data.table::fwrite(
    filtered_ccre_dt_build,
    filtered_ccre_path,
    sep = "\t"
  )
  
  message(sprintf("Total cCREs: %d", nrow(ccre_dt)))
  message(sprintf("DNase-overlapping cCREs: %d", nrow(filtered_ccre_dt_build)))
  message(sprintf("Saved overlap BED + class summary in: %s", ccre_output_dir))
}
assert_file_exists(class_summary_gz_path, "DNase-overlapping cCRE class summary TSV.GZ")
filtered_ccre_dt <- data.table::fread(filtered_ccre_path)
message(sprintf("Using filtered cCRE TSV: %s", filtered_ccre_path))

pbcpg_bed <- data.table::fread(
  cmd = sprintf("zcat '%s' | grep -v '^#'", pbcpg_bed_path),
  col.names = c("chrom", "begin", "end", "mod_score", "type", "cov",
                "est_mod_count", "est_unmod_count", "discretized_mod_score")
)[chrom == default_chrom & type == "Total"]

message(sprintf("pb-CpG BED: %d CpG sites on %s", nrow(pbcpg_bed), default_chrom))
write_gz_tsv(
  pbcpg_bed,
  file.path(ccre_output_dir, sprintf("%s_%s_pbcpg_total.tsv.gz", default_sample, default_chrom))
)

ccre_classes_base <- c("PLS", "pELS", "dELS", "CA-TF", "CA-CTCF", "CA-H3K4me3", "CA", "TF")

ccre_table_raw <- filtered_ccre_dt[
  chromosome == default_chrom & ccre_class %in% ccre_classes_base
]

if (nrow(ccre_table_raw) == 0L) {
  stop(sprintf("No DNase-overlapping cCRE rows found for chromosome: %s", default_chrom), call. = FALSE)
}

write_gz_tsv(
  ccre_table_raw,
  file.path(ccre_output_dir, sprintf("%s_%s_DNase_overlapping_cCREs.tsv.gz", default_sample, default_chrom))
)

# Synthetic Enhancers class = pELS + dELS combined
enhancer_rows <- data.table::copy(ccre_table_raw[ccre_class %in% c("pELS", "dELS")])
enhancer_rows[, ccre_class := "Enhancers"]

ccre_table_all <- data.table::rbindlist(
  list(ccre_table_raw, enhancer_rows),
  use.names = TRUE, fill = TRUE
)
ccre_table_all[, region_class := ccre_class]

write_gz_tsv(
  ccre_table_all,
  file.path(ccre_output_dir, sprintf("%s_%s_filtered_DNase_ccre_table_all.tsv.gz", default_sample, default_chrom))
)

message(sprintf("DNase-overlapping cCREs on %s: %d regions across %d classes",
                default_chrom, nrow(ccre_table_raw), length(ccre_classes_base)))

# ── Summarize all cCREs ───────────────────────────────────────────────────────
ccre_5mc_all <- data.table::rbindlist(
  lapply(seq_len(nrow(ccre_table_all)), function(i) {
    summarize_region_pbcpg_5mc(ccre_table_all[i], pbcpg_bed)
  }),
  use.names = TRUE, fill = TRUE
)

class_order <- c(ccre_classes_base, "Enhancers")
ccre_5mc_all[, region_class := factor(region_class, levels = class_order)]

class_order_dt <- data.table::data.table(
  region_class = class_order,
  class_order_idx = seq_along(class_order)
)
write_gz_tsv(
  class_order_dt,
  file.path(ccre_output_dir, sprintf("%s_%s_filtered_DNase_class_order.tsv.gz", default_sample, default_chrom))
)

write_gz_tsv(
  ccre_5mc_all,
  file.path(ccre_output_dir,
            sprintf("%s_%s_filtered_DNase_ccre_pbcpg_5mc.tsv.gz", default_sample, default_chrom))
)

message(sprintf("Regions with CpG coverage: %d / %d",
                sum(ccre_5mc_all$total_calls > 0L, na.rm = TRUE), nrow(ccre_5mc_all)))

# ── Metaprofile table only (no plotting in this chunk) ───────────────────────
meta_window_bp <- 1000L
bin_size <- 10L

ccre_wins <- data.table::copy(ccre_table_all)[, .(
  chrom = chromosome,
  region_class,
  center = as.integer((start0 + end0) / 2L),
  win_start = as.integer((start0 + end0) / 2L) - meta_window_bp,
  win_end = as.integer((start0 + end0) / 2L) + meta_window_bp
)]
ccre_wins[win_start < 0L, win_start := 0L]

cpg_for_join <- pbcpg_bed[, .(
  chrom = chrom,
  cpg_pos = begin,
  cpg_pos_end = begin,
  est_mod_count = as.numeric(est_mod_count),
  est_unmod_count = as.numeric(est_unmod_count)
)]
data.table::setkey(cpg_for_join, chrom, cpg_pos, cpg_pos_end)
data.table::setkey(ccre_wins, chrom, win_start, win_end)

overlaps <- data.table::foverlaps(
  cpg_for_join, ccre_wins,
  by.x = c("chrom", "cpg_pos", "cpg_pos_end"),
  by.y = c("chrom", "win_start", "win_end"),
  type = "within",
  nomatch = NULL
)

overlaps[, rel_pos := cpg_pos - center]
overlaps[, bin := as.integer(round(rel_pos / bin_size) * bin_size)]

meta_dt <- overlaps[, .(
  modified_calls = sum(est_mod_count),
  unmodified_calls = sum(est_unmod_count),
  total_calls = sum(est_mod_count + est_unmod_count),
  n_cpg_obs = .N
), by = .(region_class, bin)]

meta_dt[, fraction_modified := data.table::fifelse(
  total_calls > 0, modified_calls / total_calls, NA_real_
)]

all_bins <- seq(-meta_window_bp, meta_window_bp, by = bin_size)
meta_full <- data.table::CJ(
  region_class = unique(ccre_table_all$region_class),
  bin = all_bins
)
meta_full <- meta_dt[meta_full, on = .(region_class, bin)]
meta_full[, region_class := factor(region_class, levels = class_order)]

n_ccre_class <- ccre_table_all[, .N, by = .(region_class = region_class)]

# add window and bin size when write the file
# window = 1000bp, bin = 10bp

window_size <- 1000
bin_size <- 10

# Write n_ccre_class file
write_gz_tsv(
  n_ccre_class,
  file.path(ccre_output_dir, 
            sprintf("%s_%s_w%d_b%d_filtered_DNase_n_ccre_class.tsv.gz", 
                    default_sample, default_chrom, window_size, bin_size))
)

# Write metaprofile file
write_gz_tsv(
  meta_full,
  file.path(ccre_output_dir,
            sprintf("%s_%s_w%d_b%d_filtered_DNase_ccre_pbcpg_5mc_metaprofile.tsv.gz", 
                    default_sample, default_chrom, window_size, bin_size))
)

message("Saved all intermediate files for downstream plotting.")
message(sprintf("Output directory: %s", ccre_output_dir))