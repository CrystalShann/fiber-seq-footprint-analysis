source("/project/spott/cshan/fiber-seq/plot_nuc_distribution_functions.R")

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(ggplot2)
})

# ── Paths ─────────────────────────────────────────────────────────────────────
cage_tss_path <- "/project/spott/cshan/annotations/fantom5/fantom5.hg38.LCL.consensus.CAGE_peaks.withGene.bed.gz"

PEAKS_DIR <- "/project/spott/1_Shared_projects/LCL_Fiber_seq/FiberHMM/merged/combined/joint_trained_peaks/"
FP_FILES <- list(
  "10-30"   = file.path(PEAKS_DIR, "combined_all_chrs_10-30bp_fps.bed.gz"),
  "30-45"   = file.path(PEAKS_DIR, "combined_all_chrs_30-45bp_fps.bed.gz"),
  "45-60"   = file.path(PEAKS_DIR, "combined_all_chrs_45-60bp_fps.bed.gz"),
  "60-80"   = file.path(PEAKS_DIR, "combined_all_chrs_60-80bp_fps.bed.gz"),
  "140-160" = file.path(PEAKS_DIR, "combined_all_chrs_140-160bp_fps.bed.gz")
)

default_sample  <- "AL10_bc2178_19130"
window_bp       <- 2000L
bin_size        <- 10L

out_dir         <- "/project/spott/cshan/fiber-seq/results/PolII/per_read_FP"

fp_size_levels <- c("10-30", "30-45", "45-60", "60-80", "140-160")
fp_colors <- c(
  "10-30"   = "#1f77b4",
  "30-45"   = "#ff7f0e",
  "45-60"   = "#2ca02c",
  "60-80"   = "#d62728",
  "140-160" = "#9467bd"
)

# Human autosome + sex chr set used for analysis
default_chroms <- c(paste0("chr", 1:22), "chrX", "chrY")

# ── Argument parsing helpers ───────────────────────────────────────────────────
`%||%` <- function(x, y) {
  if (is.null(x) || !length(x) || is.na(x) || !nzchar(x)) y else x
}

parse_named_args <- function(args) {
  if (!length(args)) return(list())
  out <- list()
  for (arg in args) {
    if (!grepl("=", arg, fixed = TRUE)) next
    key   <- sub("=.*$", "", arg)
    value <- sub("^[^=]*=", "", arg)
    out[[key]] <- value
  }
  out
}

resolve_run_config <- function() {
  args <- parse_named_args(commandArgs(trailingOnly = TRUE))
  list(
    sample    = args[["sample"]]   %||% Sys.getenv("FIBER_SAMPLE",    default_sample),
    mode      = tolower(args[["mode"]] %||% Sys.getenv("FIBER_MODE", "combine")),
    fp_size   = args[["fp_size"]]  %||% Sys.getenv("FIBER_FP_SIZE",  ""),
    chrom     = args[["chrom"]]    %||% Sys.getenv("FIBER_CHROM",    ""),
    window_bp = as.integer(args[["window_bp"]] %||% Sys.getenv("FIBER_WINDOW_BP", as.character(window_bp))),
    bin_size  = as.integer(args[["bin_size"]]  %||% Sys.getenv("FIBER_BIN_SIZE",  as.character(bin_size)))
  )
}

# ── Caching helpers (parallel to all_chr_per_read nucleosome script) ──────────
file_is_nonempty <- function(path) file.exists(path) && isTRUE(file.info(path)$size > 0)

size_chunk_dir <- function(base_dir, sample_name, size_label) {
  file.path(base_dir, sample_name, "size_chunks", gsub("-", "_", size_label))
}

size_chunk_paths <- function(base_dir, sample_name, size_label) {
  d <- size_chunk_dir(base_dir, sample_name, size_label)
  list(
    chunk_dir = d,
    meta      = file.path(d, sprintf("per_read_%s_%s_meta.tsv.gz",    sample_name, size_label)),
    config    = file.path(d, sprintf("per_read_%s_%s_config.tsv",     sample_name, size_label)),
    partial   = file.path(d, sprintf("per_read_%s_%s_partial.rds",    sample_name, size_label))
  )
}

chunk_is_complete <- function(base_dir, sample_name, size_label) {
  p <- size_chunk_paths(base_dir, sample_name, size_label)
  file_is_nonempty(p$partial) && file_is_nonempty(p$meta)
}

# ── Load CAGE TSS ──────────────────────────────────────────────────────────────
load_cage_tss <- function(chrom_filter = NULL) {
  tss_gr <- load_tss_annotations_bed9_with_gene_id(cage_tss_path)
  if (!is.null(chrom_filter) && length(chrom_filter)) {
    tss_gr <- tss_gr[as.character(seqnames(tss_gr)) %in% chrom_filter]
  }
  tss_gr
}

# ── Read footprint BED4: chrom, start, end, RID ───────────────────────────────
# BED4 format: chrom | start (0-based) | end | read_id
# Adds a midpoint column for downstream GRanges construction.
read_fp_bed4 <- function(fp_path, chrom_filter = NULL) {
  dt <- fread(fp_path, header = FALSE, select = 1:4,
              col.names = c("chrom", "start", "end", "RID"))
  if (!is.null(chrom_filter) && length(chrom_filter)) {
    dt <- dt[chrom %in% chrom_filter]
  }
  if (!nrow(dt)) return(dt[, mid := integer()])
  dt[, mid := as.integer((start + end) / 2L)]
  dt[]
}

# ── Find footprint midpoints near TSS, tracking read identity ─────────────────
# Returns a data.table with one row per (footprint midpoint, TSS) pair within
# window_bp. Columns: RID, tss_idx, rel_pos, bin.
# The key difference from find_fp_hits_near_tss in the original script is that
# we keep RID so we can later normalize by unique reads per bin.
find_fp_hits_near_tss_per_read <- function(fp_dt, tss_gr,
                                            window_bp = 2000L,
                                            bin_size  = 10L) {
  empty <- data.table(
    RID = character(), tss_idx = integer(),
    rel_pos = integer(), bin = integer()
  )
  if (!nrow(fp_dt) || !length(tss_gr)) return(empty)

  fp_gr <- GRanges(
    seqnames = fp_dt$chrom,
    ranges   = IRanges(fp_dt$mid + 1L, fp_dt$mid + 1L),
    strand   = "*",
    RID      = fp_dt$RID
  )

  tss_windows <- promoters(tss_gr, upstream = window_bp, downstream = window_bp + 1L)
  hits <- findOverlaps(fp_gr, tss_windows)
  if (!length(hits)) return(empty)

  q_idx   <- queryHits(hits)
  s_idx   <- subjectHits(hits)
  fp_pos  <- start(fp_gr)[q_idx]
  tss_pos <- start(tss_gr)[s_idx]
  tss_str <- as.character(strand(tss_gr))[s_idx]
  rel_pos <- ifelse(tss_str == "+", fp_pos - tss_pos, tss_pos - fp_pos)

  dt <- data.table(
    RID     = fp_gr$RID[q_idx],
    tss_idx = as.integer(s_idx),
    rel_pos = as.integer(rel_pos),
    bin     = as.integer(round(rel_pos / bin_size) * bin_size)
  )
  dt[bin >= -window_bp & bin <= window_bp]
}

# ── Build per-read meta profile for one footprint size class ──────────────────
# Returns a data.table with columns:
#   bin, N (raw count), n_reads (unique RIDs), n_tss (unique TSS indices),
#   N_per_read (N / n_reads per bin), N_per_tss (N / total_tss),
#   fp_size, sample, chrom_scope.
# N_per_read: coverage-normalized signal — removes bias from loci with many
#   overlapping reads. N_per_tss: matches the normalisation in the original
#   (non-per-read) script so results are directly comparable.
make_fp_meta_per_read <- function(hit_dt, n_tss_total, size_label,
                                   sample_name, scope_label,
                                   window_bp = 2000L, bin_size = 10L) {
  template <- data.table(bin = seq(-window_bp, window_bp, by = bin_size))

  counts <- if (nrow(hit_dt)) {
    hit_dt[, .(
      N       = .N,
      n_reads = uniqueN(RID),
      n_tss   = uniqueN(tss_idx)
    ), by = bin]
  } else {
    data.table(bin = integer(), N = integer(), n_reads = integer(), n_tss = integer())
  }

  meta <- merge(template, counts, by = "bin", all.x = TRUE, sort = TRUE)
  meta[is.na(N), `:=`(N = 0L, n_reads = 0L, n_tss = 0L)]

  # per-read normalization: divides by unique reads contributing to each bin
  meta[, N_per_read := fifelse(n_reads > 0L, N / n_reads, 0)]

  # per-TSS normalization: matches original script for direct comparison
  meta[, N_per_tss := if (n_tss_total > 0L) N / n_tss_total else 0]

  meta[, `:=`(fp_size = size_label, sample = sample_name, chrom_scope = scope_label)]
  meta[]
}

# ── Set factor levels ──────────────────────────────────────────────────────────
set_profile_factors <- function(dt) {
  if (!nrow(dt)) return(dt)
  if ("fp_size" %in% names(dt)) {
    dt[, fp_size := factor(fp_size, levels = fp_size_levels)]
  }
  dt
}

# ── Build plots (same layout as footprints_TSS_profiles_with_fiber_features.r) ─
# raw_col: column for raw-count panel (N)
# norm_col: column for normalized panel (N_per_read or N_per_tss)
build_scope_plots <- function(meta, sample_name, scope_label,
                               norm_col  = "N_per_read",
                               window_bp = 2000L) {
  scope_subtitle <- sprintf("Sample: %s | %s (per-read)",
                            sample_name,
                            if (identical(scope_label, "all_chr")) "All chromosomes"
                            else paste("Chromosome:", scope_label))

  x_scale <- scale_x_continuous(
    breaks = seq(-window_bp, window_bp, 500),
    labels = function(x) paste0(x / 1000, " kb")
  )

  # 1. Raw counts, faceted by fp_size
  raw_plot <- ggplot(meta, aes(x = bin, y = N, color = fp_size)) +
    geom_line(linewidth = 0.7) +
    geom_smooth(method = "loess", span = 0.06, se = FALSE, linewidth = 0.9) +
    facet_wrap(~fp_size, scales = "free_y", ncol = 1) +
    scale_color_manual(values = fp_colors, drop = FALSE) +
    x_scale +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
    labs(
      x        = "Distance to TSS",
      y        = "Footprint midpoints per 10 bp bin",
      color    = "Footprint size",
      title    = "Footprint size profiles around CAGE TSS: raw counts (per-read)",
      subtitle = scope_subtitle
    ) +
    theme_bw(base_size = 12)

  # 2. Per-read normalized, faceted by fp_size
  norm_label <- if (norm_col == "N_per_read") {
    "Mean footprint midpoints per contributing read (10 bp bins)"
  } else {
    "Average footprint midpoints per TSS (10 bp bins)"
  }
  norm_title <- if (norm_col == "N_per_read") {
    "Footprint size profiles around CAGE TSS: per-read normalized"
  } else {
    "Footprint size profiles around CAGE TSS: per-TSS normalized"
  }

  norm_plot <- ggplot(meta, aes_string(x = "bin", y = norm_col, color = "fp_size")) +
    geom_line(linewidth = 0.7) +
    geom_smooth(method = "loess", span = 0.06, se = FALSE, linewidth = 0.9) +
    facet_wrap(~fp_size, scales = "free_y", ncol = 1) +
    scale_color_manual(values = fp_colors, drop = FALSE) +
    x_scale +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
    labs(
      x        = "Distance to TSS",
      y        = norm_label,
      color    = "Footprint size",
      title    = norm_title,
      subtitle = scope_subtitle
    ) +
    theme_bw(base_size = 12)

  # 3. Overlay – all size classes on one panel (normalized)
  overlay_plot <- ggplot(meta, aes_string(x = "bin", y = norm_col, color = "fp_size")) +
    geom_line(linewidth = 0.8) +
    scale_color_manual(values = fp_colors, drop = FALSE) +
    x_scale +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
    labs(
      x        = "Distance to TSS",
      y        = norm_label,
      color    = "Footprint size",
      title    = "Footprint size profiles around CAGE TSS: all sizes overlaid (per-read)",
      subtitle = scope_subtitle
    ) +
    theme_bw(base_size = 12)

  # 4. Z-score overlay – highlights relative shape differences
  meta_z <- copy(meta)
  norm_vals <- meta_z[[norm_col]]
  meta_z[, z := {
    v <- get(norm_col)
    s <- sd(v)
    if (is.na(s) || s == 0) rep(0, .N) else as.numeric(scale(v))
  }, by = fp_size]

  zscore_plot <- ggplot(meta_z, aes(x = bin, y = z, color = fp_size)) +
    geom_line(linewidth = 0.8) +
    scale_color_manual(values = fp_colors, drop = FALSE) +
    x_scale +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
    labs(
      x        = "Distance to TSS",
      y        = "Z-scored signal",
      color    = "Footprint size",
      title    = "Footprint size profiles around CAGE TSS: Z-score overlay (per-read)",
      subtitle = "Z-scores computed within each footprint size class"
    ) +
    theme_bw(base_size = 12)

  list(
    raw_plot     = raw_plot,
    norm_plot    = norm_plot,
    overlay_plot = overlay_plot,
    zscore_plot  = zscore_plot
  )
}

# ── Write plots and tables for a scope ────────────────────────────────────────
write_scope_outputs <- function(sample_name, scope_label, meta,
                                 analysis_config, plots) {
  plot_dir  <- file.path(out_dir, "plots")
  table_dir <- file.path(out_dir, "tables")
  dir.create(plot_dir,  recursive = TRUE, showWarnings = FALSE)
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

  print(plots$raw_plot)
  print(plots$norm_plot)
  print(plots$overlay_plot)
  print(plots$zscore_plot)

  tag <- sprintf("%s_%s", sample_name, scope_label)
  ggsave(file.path(plot_dir, sprintf("%s_fp_per_read_profiles_raw.pdf",     tag)),
         plot = plots$raw_plot,     width = 8, height = 14)
  ggsave(file.path(plot_dir, sprintf("%s_fp_per_read_profiles_norm.pdf",    tag)),
         plot = plots$norm_plot,    width = 8, height = 14)
  ggsave(file.path(plot_dir, sprintf("%s_fp_per_read_profiles_overlay.pdf", tag)),
         plot = plots$overlay_plot, width = 8, height = 5)
  ggsave(file.path(plot_dir, sprintf("%s_fp_per_read_profiles_zscore.pdf",  tag)),
         plot = plots$zscore_plot,  width = 8, height = 5)

  fwrite(meta,            file.path(table_dir, sprintf("%s_fp_per_read_meta.tsv.gz",   tag)), sep = "\t")
  fwrite(analysis_config, file.path(table_dir, sprintf("%s_fp_per_read_config.tsv",    tag)), sep = "\t")

  message(sprintf("Outputs written to %s", out_dir))
}

# ── Per-size-class chunk ───────────────────────────────────────────────────────
# Processes one footprint size class: loads BED4, finds hits near TSS per read,
# builds meta profile, saves a partial RDS for later combine step.
run_size_chunk <- function(config, tss_gr, size_label) {
  if (!size_label %in% names(FP_FILES)) {
    stop("Unknown fp_size: ", size_label, ". Valid options: ", paste(names(FP_FILES), collapse = ", "))
  }

  out_base <- file.path(out_dir, "chunks")
  paths    <- size_chunk_paths(out_base, config$sample, size_label)

  if (chunk_is_complete(out_base, config$sample, size_label)) {
    message("  Reusing cached chunk for size class: ", size_label)
    return(invisible(readRDS(paths$partial)))
  }

  dir.create(paths$chunk_dir, recursive = TRUE, showWarnings = FALSE)

  fp_path <- FP_FILES[[size_label]]
  message(sprintf("  Loading footprint BED4 [%s]: %s", size_label, basename(fp_path)))

  # Optional chrom filter: restrict to requested chromosomes to save memory
  chrom_filter <- if (nzchar(config$chrom)) {
    trimws(unlist(strsplit(config$chrom, ",", fixed = TRUE)))
  } else {
    default_chroms
  }

  fp_dt <- read_fp_bed4(fp_path, chrom_filter = chrom_filter)
  message(sprintf("    %d footprints loaded", nrow(fp_dt)))

  hit_dt <- find_fp_hits_near_tss_per_read(
    fp_dt, tss_gr,
    window_bp = config$window_bp,
    bin_size  = config$bin_size
  )
  message(sprintf("    %d footprint-TSS hit pairs within +/-%d bp", nrow(hit_dt), config$window_bp))
  rm(fp_dt); gc(verbose = FALSE)

  meta <- make_fp_meta_per_read(
    hit_dt,
    n_tss_total = length(tss_gr),
    size_label  = size_label,
    sample_name = config$sample,
    scope_label = "all_chr",
    window_bp   = config$window_bp,
    bin_size    = config$bin_size
  )

  fwrite(meta,
         paths$meta,
         sep = "\t")
  fwrite(
    data.table(
      sample    = config$sample,
      fp_size   = size_label,
      chrom     = config$chrom %||% "all_chr",
      window_bp = config$window_bp,
      bin_size  = config$bin_size,
      n_tss     = length(tss_gr),
      n_fp_hits = nrow(hit_dt)
    ),
    paths$config,
    sep = "\t"
  )

  obj <- list(sample = config$sample, fp_size = size_label, meta = meta)
  saveRDS(obj, paths$partial)

  message(sprintf("  Chunk saved: %s", paths$partial))
  invisible(obj)
}

# ── Combine all size-class chunks into final plots ─────────────────────────────
combine_size_chunks <- function(config) {
  out_base <- file.path(out_dir, "chunks")

  missing_sizes <- fp_size_levels[!vapply(
    fp_size_levels,
    function(s) chunk_is_complete(out_base, config$sample, s),
    logical(1)
  )]
  if (length(missing_sizes)) {
    stop(
      "Missing chunk files. Run mode=size_chunk for these fp_size classes first:\n",
      paste(missing_sizes, collapse = "\n")
    )
  }

  metas <- rbindlist(
    lapply(fp_size_levels, function(s) {
      readRDS(size_chunk_paths(out_base, config$sample, s)$partial)$meta
    }),
    use.names = TRUE, fill = TRUE
  )
  metas <- set_profile_factors(metas)

  # Use the n_tss from the first chunk's config as the global TSS count
  first_config_path <- size_chunk_paths(out_base, config$sample, fp_size_levels[1])$config
  n_tss_global <- if (file.exists(first_config_path)) {
    fread(first_config_path)$n_tss[1]
  } else {
    NA_integer_
  }

  scope_label <- "all_chr"
  plots <- build_scope_plots(
    meta        = metas,
    sample_name = config$sample,
    scope_label = scope_label,
    norm_col    = "N_per_read",
    window_bp   = config$window_bp
  )

  analysis_config <- data.table(
    sample      = config$sample,
    chrom_scope = scope_label,
    n_tss       = n_tss_global,
    tss_source  = "CAGE",
    window_bp   = config$window_bp,
    bin_size    = config$bin_size,
    fp_sizes    = paste(fp_size_levels, collapse = ","),
    norm_method = "N_per_read"
  )

  write_scope_outputs(
    sample_name     = config$sample,
    scope_label     = scope_label,
    meta            = metas,
    analysis_config = analysis_config,
    plots           = plots
  )

  invisible(metas)
}

# ── Entry point ───────────────────────────────────────────────────────────────
run_footprint_TSS_profiles_per_read <- function() {
  config <- resolve_run_config()
  message(
    "mode=", config$mode,
    " | sample=", config$sample,
    " | window_bp=", config$window_bp,
    " | bin_size=", config$bin_size
  )

  # mode=size_chunk: process one fp_size class, save partial RDS
  if (config$mode == "size_chunk") {
    if (!nzchar(config$fp_size)) {
      stop("mode=size_chunk requires fp_size=<size>, e.g. fp_size=10-30")
    }
    tss_gr <- load_cage_tss()
    message("CAGE TSS loaded: ", length(tss_gr), " sites")
    return(invisible(run_size_chunk(config, tss_gr, config$fp_size)))
  }

  # mode=combine: merge all existing size-class chunks and plot
  if (config$mode == "combine") {
    return(invisible(combine_size_chunks(config)))
  }

  # mode=full: run all size-class chunks then combine
  if (config$mode == "full") {
    tss_gr <- load_cage_tss()
    message("CAGE TSS loaded: ", length(tss_gr), " sites")
    for (size_label in fp_size_levels) {
      message("Processing size class: ", size_label)
      run_size_chunk(config, tss_gr, size_label)
      gc(verbose = FALSE)
    }
    return(invisible(combine_size_chunks(config)))
  }

  stop("Unsupported mode: ", config$mode,
       ". Use mode=size_chunk (requires fp_size=<size>), mode=combine, or mode=full.")
}

run_footprint_TSS_profiles_per_read()

# ── Usage examples ──────────────────────────────────────────────────────────────
# Process one size class at a time (submit as SLURM array, one job per size):
#   Rscript /project/spott/cshan/fiber-seq/footprints_TSS_profiles_per_read.r \
#     sample=AL10_bc2178_19130 mode=size_chunk fp_size=10-30
#   Rscript /project/spott/cshan/fiber-seq/footprints_TSS_profiles_per_read.r \
#     sample=AL10_bc2178_19130 mode=size_chunk fp_size=30-45
#   Rscript /project/spott/cshan/fiber-seq/footprints_TSS_profiles_per_read.r \
#     sample=AL10_bc2178_19130 mode=size_chunk fp_size=45-60
#   Rscript /project/spott/cshan/fiber-seq/footprints_TSS_profiles_per_read.r \
#     sample=AL10_bc2178_19130 mode=size_chunk fp_size=60-80
#   Rscript /project/spott/cshan/fiber-seq/footprints_TSS_profiles_per_read.r \
#     sample=AL10_bc2178_19130 mode=size_chunk fp_size=140-160
#
# After all size_chunk jobs complete, combine and plot:
#   Rscript /project/spott/cshan/fiber-seq/footprints_TSS_profiles_per_read.r \
#     sample=AL10_bc2178_19130 mode=combine
#
# Or run everything in one job (memory-intensive):
#   Rscript /project/spott/cshan/fiber-seq/footprints_TSS_profiles_per_read.r \
#     sample=AL10_bc2178_19130 mode=full
#
# To restrict to specific chromosomes during size_chunk (e.g. for testing):
#   Rscript /project/spott/cshan/fiber-seq/footprints_TSS_profiles_per_read.r \
#     sample=AL10_bc2178_19130 mode=size_chunk fp_size=10-30 chrom=chr1
#
#SBATCH --job-name=fp_per_read_tss
#SBATCH --time=24:00:00
#SBATCH --mem=120G
#SBATCH --cpus-per-task=3
#SBATCH --account=pi-spott
#SBATCH --partition=bigmem
#SBATCH --error=/project/spott/cshan/fiber-seq/results/logs/fp_per_read_tss_%A_%a.err
