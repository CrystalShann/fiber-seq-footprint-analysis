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
default_sample  <- "AL10_bc2178_19130"
default_chrom   <- "chr1"

results_root_dir <- "/project/spott/cshan/fiber-seq/results/promoter_summary/CpG_methylation_cCRE"
ccre_output_dir  <- file.path(results_root_dir, "ccre")

pbcpg_bed_path <- file.path(
  "/project/spott/cshan/fiber-seq/results/fire_CpG",
  default_sample,
  sprintf("%s_CPG.combined.bed.gz", default_sample)
)
ccre_path <- "/project/spott/cshan/annotations/GRCh38-cCREs.bed"

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
      end   <= region_row$end0
  ]
  modified_calls   <- as.integer(round(sum(hits$est_mod_count)))
  unmodified_calls <- as.integer(round(sum(hits$est_unmod_count)))
  total_calls      <- modified_calls + unmodified_calls
  data.table::data.table(
    region_class            = region_row$region_class,
    chromosome              = region_row$chromosome,
    start                   = region_row$start0 + 1L,
    end                     = region_row$end0,
    region_width_bp         = region_row$end0 - region_row$start0,
    n_cpg_sites             = nrow(hits),
    modified_calls          = modified_calls,
    unmodified_calls        = unmodified_calls,
    total_calls             = total_calls,
    fraction_modified       = data.table::fifelse(
      total_calls > 0L, modified_calls / total_calls, NA_real_
    )
  )
}

# ── Load data ─────────────────────────────────────────────────────────────────
assert_file_exists(pbcpg_bed_path, "pb-CpG-tools combined BED")
assert_file_exists(ccre_path,      "cCRE BED")

pbcpg_bed <- data.table::fread(
  cmd = sprintf("zcat '%s' | grep -v '^#'", pbcpg_bed_path),
  col.names = c("chrom", "begin", "end", "mod_score", "type", "cov",
                "est_mod_count", "est_unmod_count", "discretized_mod_score")
)[chrom == default_chrom & type == "Total"]

message(sprintf("pb-CpG BED: %d CpG sites on %s", nrow(pbcpg_bed), default_chrom))

ccre_classes_base <- c("PLS", "pELS", "dELS", "CA-TF", "CA-CTCF", "CA-H3K4me3", "CA", "TF")

ccre_table_raw <- data.table::fread(
  ccre_path, header = FALSE,
  col.names = c("chromosome", "start0", "end0", "ccre_id", "ccre_element_id", "ccre_class")
)[chromosome == default_chrom & ccre_class %in% ccre_classes_base]

# Synthetic Enhancers class = pELS + dELS combined
enhancer_rows <- data.table::copy(ccre_table_raw[ccre_class %in% c("pELS", "dELS")])
enhancer_rows[, ccre_class := "Enhancers"]

ccre_table_all <- data.table::rbindlist(
  list(ccre_table_raw, enhancer_rows),
  use.names = TRUE, fill = TRUE
)
ccre_table_all[, region_class := ccre_class]

message(sprintf("cCREs loaded: %d regions across %d classes on %s",
                nrow(ccre_table_raw), length(ccre_classes_base), default_chrom))

# ── Summarize all cCREs ───────────────────────────────────────────────────────
ccre_5mc_all <- data.table::rbindlist(
  lapply(seq_len(nrow(ccre_table_all)), function(i) {
    summarize_region_pbcpg_5mc(ccre_table_all[i], pbcpg_bed)
  }),
  use.names = TRUE, fill = TRUE
)

# Order classes: base classes in defined order, then synthetic Enhancers
class_order <- c(ccre_classes_base, "Enhancers")
ccre_5mc_all[, region_class := factor(region_class, levels = class_order)]

message(sprintf("Regions with CpG coverage: %d / %d",
                sum(ccre_5mc_all$total_calls > 0L, na.rm = TRUE), nrow(ccre_5mc_all)))

# ── Plot ──────────────────────────────────────────────────────────────────────
n_per_class <- ccre_5mc_all[!is.na(fraction_modified), .N, by = region_class]
x_labels    <- setNames(
  sprintf("%s\n(n=%d)", n_per_class$region_class, n_per_class$N),
  n_per_class$region_class
)

all_ccre_5mc_plot <- ggplot2::ggplot(
  ccre_5mc_all[!is.na(fraction_modified)],
  ggplot2::aes(x = region_class, y = fraction_modified)
) +
  ggplot2::geom_boxplot(fill = "#8C510A", alpha = 0.45, outlier.shape = NA) +
  ggplot2::geom_jitter(width = 0.2, size = 0.6, alpha = 0.25, color = "#8C510A") +
  ggplot2::scale_x_discrete(labels = x_labels) +
  ggplot2::scale_y_continuous(limits = c(0, 1), expand = ggplot2::expansion(mult = 0.02)) +
  ggplot2::labs(
    title = sprintf("CpG-only 5mC fraction by cCRE class — all %s cCREs (%s)",
                    default_chrom, default_sample),
    x     = "cCRE class",
    y     = "CpG methylation fraction (pb-CpG-tools)"
  ) +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(size = 10))

# ── Save ──────────────────────────────────────────────────────────────────────
write_gz_tsv(
  ccre_5mc_all,
  file.path(ccre_output_dir,
            sprintf("%s_%s_all_ccre_pbcpg_5mc.tsv.gz", default_sample, default_chrom))
)

ggplot2::ggsave(
  file.path(ccre_output_dir,
            sprintf("%s_%s_all_ccre_pbcpg_5mc.pdf", default_sample, default_chrom)),
  all_ccre_5mc_plot, width = 14, height = 6
)



# fwrite(x = pbcpg_bed, 
#        "/project/spott/cshan/fiber-seq/results/promoter_summary/CpG_methylation_cCRE/pbcpg_bed_chr1_AL10.csv")
# fwrite(x = ccre_table_all, 
#        "/project/spott/cshan/fiber-seq/results/promoter_summary/CpG_methylation_cCRE/ccre_table_all_chr1_AL10.csv")





## All-cCRE CpG-only 5mC metaprofile (aggregate signal around cCRE centers)

# Average CpG methylation fraction binned by distance from each cCRE center.
# Uses all chr1 cCREs from the previous chunk (`pbcpg_bed`, `ccre_table_all`, `class_order`).



meta_window_bp <- 300L   # flanking bp each side of cCRE center
bin_size       <- 1L     # bp per bin

# ── Build windows centred on each cCRE ───────────────────────────────────────
ccre_wins <- data.table::copy(ccre_table_all)[, .(
  chrom      = chromosome,
  region_class,
  center     = as.integer((start0 + end0) / 2L),
  win_start  = as.integer((start0 + end0) / 2L) - meta_window_bp,
  win_end    = as.integer((start0 + end0) / 2L) + meta_window_bp
)]
ccre_wins[win_start < 0L, win_start := 0L]   # clamp to chromosome start

# ── Prepare CpG table for interval overlap ────────────────────────────────────
# pbcpg_bed already filtered to default_chrom & type == "Total"
cpg_for_join <- pbcpg_bed[, .(
  chrom           = chrom,
  cpg_pos         = begin,           # 0-based; treat as point interval
  cpg_pos_end     = begin,
  est_mod_count   = as.numeric(est_mod_count),
  est_unmod_count = as.numeric(est_unmod_count)
)]
data.table::setkey(cpg_for_join, chrom, cpg_pos, cpg_pos_end)
data.table::setkey(ccre_wins,    chrom, win_start, win_end)

# foverlaps: for each CpG point, find all cCRE windows it falls in
overlaps <- data.table::foverlaps(
  cpg_for_join, ccre_wins,
  by.x    = c("chrom", "cpg_pos", "cpg_pos_end"),
  by.y    = c("chrom", "win_start", "win_end"),
  type    = "within",
  nomatch = NULL
)

overlaps[, rel_pos := cpg_pos - center]
overlaps[, bin     := as.integer(round(rel_pos / bin_size) * bin_size)]

# ── Aggregate per (class, bin) ────────────────────────────────────────────────
meta_dt <- overlaps[, .(
  modified_calls   = sum(est_mod_count),
  unmodified_calls = sum(est_unmod_count),
  total_calls      = sum(est_mod_count + est_unmod_count),
  n_cpg_obs        = .N
), by = .(region_class, bin)]

meta_dt[, fraction_modified := data.table::fifelse(
  total_calls > 0, modified_calls / total_calls, NA_real_
)]

# Fill in missing bins so lines are continuous
all_bins  <- seq(-meta_window_bp, meta_window_bp, by = bin_size)
meta_full <- data.table::CJ(
  region_class = unique(ccre_table_all$region_class),
  bin          = all_bins
)
meta_full <- meta_dt[meta_full, on = .(region_class, bin)]
meta_full[, region_class := factor(region_class, levels = class_order)]

# Count of cCREs per class (for facet labels)
n_ccre_class <- ccre_table_all[, .N, by = .(region_class = region_class)]
facet_labels <- setNames(
  sprintf("%s\n(n=%d cCREs)", n_ccre_class$region_class, n_ccre_class$N),
  n_ccre_class$region_class
)

# ── Plot ──────────────────────────────────────────────────────────────────────
meta_plot <- ggplot2::ggplot(
  meta_full[!is.na(fraction_modified)],
  ggplot2::aes(x = bin, y = fraction_modified)
) +
  ggplot2::geom_line(color = "#8C510A", linewidth = 0.4, na.rm = TRUE) +
  ggplot2::geom_smooth(
    method    = "loess",
    span      = 0.15,
    se        = FALSE,
    color     = "#8C510A",
    linewidth = 0.9,
    na.rm     = TRUE
  ) +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
                      color = "grey40", linewidth = 0.4) +
  ggplot2::scale_x_continuous(
    breaks = seq(-meta_window_bp, meta_window_bp, 500),
    labels = function(x) paste0(x / 1000, " kb")
  ) +
  ggplot2::scale_y_continuous(
    limits = c(0, 1),
    expand = ggplot2::expansion(mult = 0.02)
  ) +
  ggplot2::facet_wrap(
    ~region_class,
    nrow    = 3,
    labeller = ggplot2::as_labeller(facet_labels)
  ) +
  ggplot2::labs(
    title = sprintf(
      "CpG-only 5mC metaprofile around cCRE centers — all %s cCREs (%s)",
      default_chrom, default_sample
    ),
    x = "Distance from cCRE center",
    y = "CpG methylation fraction (pb-CpG-tools)"
  ) +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(
    strip.text   = ggplot2::element_text(size = 9),
    axis.text.x  = ggplot2::element_text(size = 8)
  )

# ── Save ──────────────────────────────────────────────────────────────────────
write_gz_tsv(
  meta_full,
  file.path(ccre_output_dir,
            sprintf("%s_%s_all_ccre_pbcpg_5mc_metaprofile.tsv.gz", default_sample, default_chrom))
)

ggplot2::ggsave(
  file.path(ccre_output_dir,
            sprintf("%s_%s_all_ccre_pbcpg_5mc_metaprofile.pdf", default_sample, default_chrom)),
  meta_plot, width = 16, height = 10
)
