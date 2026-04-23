library(data.table)
library(ggplot2)
library(patchwork)
library(scales)
library(zoo)

# ── Variables ─────────────────────────────────────────────────────────────────
default_sample <- "AL10_bc2178_19130"
meta_window_bp <- 1000L
bin_size       <- 10L
smooth_window  <- 5L

all_chroms <- paste0("chr", c(1:22, "X", "Y"))

output_dir        <- "/project/spott/cshan/fiber-seq/results/cCRE_summary/silencer_summary/CpG_silencers_v2"
silencer_path_raw <- "/project/spott/cshan/annotations/human_silencers/STARR-Silencers.Stringent.bed"
pbcpg_bed_path    <- file.path(
  "/project/spott/cshan/fiber-seq/results/fire_CpG",
  default_sample,
  sprintf("%s_CPG.combined.bed.gz", default_sample)
)

# ── Load raw silencer BED (all chroms) ────────────────────────────────────────
message("Loading silencer BED...")
silencer_dt <- data.table::fread(
  silencer_path_raw,
  header    = FALSE,
  col.names = c("chromosome", "start0", "end0", "id1", "id2", "ccre_class")
)
silencer_dt <- silencer_dt[chromosome %in% all_chroms]
silencer_dt[, silencer_id := .I]

n_per_class <- silencer_dt[, .N, by = ccre_class]
ccre_classes <- sort(unique(silencer_dt$ccre_class))
message(sprintf("Total silencers (all chroms): %d", nrow(silencer_dt)))
message(sprintf("Classes: %s",
                paste(sprintf("%s(n=%d)", n_per_class$ccre_class, n_per_class$N), collapse = ", ")))

# Save silencer summary
data.table::fwrite(
  silencer_dt,
  file.path(output_dir, sprintf("%s_all_chroms_silencers_with_class.tsv.gz", default_sample)),
  sep = "\t"
)

# ── Load raw pbcpg BED (all chroms) ───────────────────────────────────────────
message("Loading pbcpg BED (all chroms, this may take a moment)...")
pbcpg_bed <- data.table::fread(
  cmd       = sprintf("zcat '%s' | grep -v '^#'", pbcpg_bed_path),
  col.names = c("chrom", "begin", "end", "mod_score", "type", "cov",
                "est_mod_count", "est_unmod_count", "discretized_mod_score")
)[chrom %in% all_chroms & type == "Total"]

message(sprintf("pb-CpG BED: %d CpG sites across all chromosomes", nrow(pbcpg_bed)))

# ── Build metaprofile windows ─────────────────────────────────────────────────
message("Computing metaprofile windows...")
silencer_wins <- data.table::copy(silencer_dt)[, .(
  chrom       = chromosome,
  silencer_id,
  ccre_class,
  center      = as.integer((start0 + end0) / 2L),
  win_start   = as.integer((start0 + end0) / 2L) - meta_window_bp,
  win_end     = as.integer((start0 + end0) / 2L) + meta_window_bp
)]
silencer_wins[win_start < 0L, win_start := 0L]

cpg_for_join <- pbcpg_bed[, .(
  chrom           = chrom,
  cpg_pos         = begin,
  cpg_pos_end     = begin,
  est_mod_count   = as.numeric(est_mod_count),
  est_unmod_count = as.numeric(est_unmod_count)
)]
data.table::setkey(cpg_for_join,  chrom, cpg_pos, cpg_pos_end)
data.table::setkey(silencer_wins, chrom, win_start, win_end)

message("Running foverlaps (all chroms x all silencers)...")
overlaps <- data.table::foverlaps(
  cpg_for_join, silencer_wins,
  by.x    = c("chrom", "cpg_pos", "cpg_pos_end"),
  by.y    = c("chrom", "win_start", "win_end"),
  type    = "within",
  nomatch = NULL
)

overlaps[, rel_pos := cpg_pos - center]
overlaps[, bin     := as.integer(round(rel_pos / bin_size) * bin_size)]

# ── Aggregate per-silencer methylation for boxplot ───────────────────────────
message("Aggregating per-silencer methylation fractions...")
silencer_5mc_all <- overlaps[, .(
  modified_calls   = sum(est_mod_count),
  unmodified_calls = sum(est_unmod_count),
  total_calls      = sum(est_mod_count + est_unmod_count),
  n_cpg_sites      = uniqueN(cpg_pos)
), by = .(silencer_id, ccre_class, chromosome, start0, end0)]

silencer_5mc_all <- silencer_dt[, .(
  silencer_id,
  chromosome,
  start0,
  end0,
  ccre_class
)][silencer_5mc_all, on = .(silencer_id, ccre_class, chromosome, start0, end0)]

silencer_5mc_all[is.na(modified_calls), `:=`(
  modified_calls   = 0,
  unmodified_calls = 0,
  total_calls      = 0,
  n_cpg_sites      = 0
)]
silencer_5mc_all[, region_width_bp := end0 - start0]
silencer_5mc_all[, fraction_modified := data.table::fifelse(
  total_calls > 0, modified_calls / total_calls, NA_real_
)]
silencer_5mc_all[, ccre_class := factor(ccre_class, levels = ccre_classes)]

data.table::fwrite(
  silencer_5mc_all,
  file.path(output_dir, sprintf("%s_all_chroms_silencer_pbcpg_5mc.tsv.gz", default_sample)),
  sep = "\t"
)

# ── Aggregate metaprofile: all silencers combined ─────────────────────────────
message("Aggregating metaprofile (all silencers)...")
meta_all <- overlaps[, .(
  modified_calls   = sum(est_mod_count),
  unmodified_calls = sum(est_unmod_count),
  total_calls      = sum(est_mod_count + est_unmod_count),
  n_cpg_obs        = .N
), by = .(bin)]

meta_all[, fraction_modified := data.table::fifelse(
  total_calls > 0, modified_calls / total_calls, NA_real_
)]

all_bins  <- seq(-meta_window_bp, meta_window_bp, by = bin_size)
meta_all_full <- data.table::data.table(bin = all_bins)
meta_all_full <- meta_all[meta_all_full, on = .(bin)]

data.table::fwrite(
  meta_all_full,
  file.path(output_dir,
            sprintf("%s_all_chroms_w%d_b%d_silencer_pbcpg_5mc_metaprofile.tsv.gz",
                    default_sample, meta_window_bp, bin_size)),
  sep = "\t"
)

# ── Aggregate metaprofile: per cCRE class ─────────────────────────────────────
message("Aggregating metaprofile by cCRE class...")
meta_class <- overlaps[, .(
  modified_calls   = sum(est_mod_count),
  unmodified_calls = sum(est_unmod_count),
  total_calls      = sum(est_mod_count + est_unmod_count),
  n_cpg_obs        = .N
), by = .(ccre_class, bin)]

meta_class[, fraction_modified := data.table::fifelse(
  total_calls > 0, modified_calls / total_calls, NA_real_
)]

meta_class_full <- data.table::CJ(
  ccre_class = ccre_classes,
  bin        = all_bins
)
meta_class_full <- meta_class[meta_class_full, on = .(ccre_class, bin)]

data.table::fwrite(
  meta_class_full,
  file.path(output_dir,
            sprintf("%s_all_chroms_w%d_b%d_silencer_pbcpg_5mc_metaprofile_by_class.tsv.gz",
                    default_sample, meta_window_bp, bin_size)),
  sep = "\t"
)

# ── Smooth ────────────────────────────────────────────────────────────────────
meta_all_full[, fraction_modified_smooth := zoo::rollmean(
  fraction_modified, k = smooth_window, fill = NA, align = "center"
)]

meta_class_full[, fraction_modified_smooth := zoo::rollmean(
  fraction_modified, k = smooth_window, fill = NA, align = "center"
), by = ccre_class]

n_lookup <- setNames(n_per_class$N, n_per_class$ccre_class)
meta_class_full[, class_label := paste0(
  ccre_class, "\n(n=", n_lookup[ccre_class], ")"
)]

# ── Colors ────────────────────────────────────────────────────────────────────
class_colors <- c(
  "PLS"        = "#d73027",
  "pELS"       = "#fc8d59",
  "dELS"       = "#fee090",
  "CA-TF"      = "#91bfdb",
  "CA-CTCF"    = "#4575b4",
  "CA-H3K4me3" = "#1a9850",
  "CA"         = "#a6d96a",
  "TF"         = "#762a83"
)
missing_cols <- setdiff(ccre_classes, names(class_colors))
if (length(missing_cols) > 0) {
  extra <- scales::hue_pal()(length(missing_cols))
  names(extra) <- missing_cols
  class_colors <- c(class_colors, extra)
}

# ── Plot 0: boxplot of per-silencer methylation by cCRE class ────────────────
message("Plotting per-silencer methylation boxplot...")
n_box <- silencer_5mc_all[!is.na(fraction_modified), .N, by = ccre_class]
x_labels <- setNames(
  sprintf("%s\n(n=%d)", n_box$ccre_class, n_box$N),
  n_box$ccre_class
)

p_box <- ggplot2::ggplot(
  silencer_5mc_all[!is.na(fraction_modified)],
  ggplot2::aes(x = ccre_class, y = fraction_modified)
) +
  ggplot2::geom_boxplot(fill = "#8C510A", alpha = 0.45, outlier.shape = NA) +
  ggplot2::geom_jitter(
    width = 0.2, size = 0.6, alpha = 0.25, color = "#8C510A", na.rm = TRUE
  ) +
  ggplot2::scale_x_discrete(labels = x_labels, drop = FALSE) +
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = 0.02)) +
  ggplot2::labs(
    title = sprintf(
      "CpG-only 5mC fraction by silencer cCRE class - all chromosomes (%s)",
      default_sample
    ),
    x = "cCRE class",
    y = "CpG methylation fraction (pb-CpG-tools)"
  ) +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(size = 10))

ggplot2::ggsave(
  file.path(output_dir,
            sprintf("%s_all_chroms_silencer_pbcpg_5mc_boxplot.pdf", default_sample)),
  p_box, width = 14, height = 6, device = "pdf"
)
message("Saved boxplot.")

# ── Plot 1: all silencers combined ────────────────────────────────────────────
message("Plotting all silencers combined...")
n_total   <- nrow(silencer_dt)
plot_all  <- meta_all_full[!is.na(fraction_modified_smooth)]

p_all_meta <- ggplot2::ggplot(plot_all, ggplot2::aes(x = bin, y = fraction_modified_smooth)) +
  ggplot2::geom_line(color = "#2166ac", linewidth = 0.8) +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.5) +
  ggplot2::scale_x_continuous(
    breaks = seq(-meta_window_bp, meta_window_bp, by = 250),
    labels = function(x) paste0(x, " bp")
  ) +
  ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  ggplot2::labs(
    title    = sprintf("5mCpG Metaprofile at STARR-Silencers — All Classes (n = %d)", n_total),
    subtitle = sprintf("Sample: %s | All chromosomes | Window: ±%d bp | Bin: %d bp",
                       default_sample, meta_window_bp, bin_size),
    x = "Distance from silencer center (bp)",
    y = "Fraction 5mCpG"
  ) +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::theme(
    plot.title       = ggplot2::element_text(face = "bold", size = 13),
    plot.subtitle    = ggplot2::element_text(size = 10, color = "grey40"),
    panel.grid.minor = ggplot2::element_blank(),
    axis.text.x      = ggplot2::element_text(angle = 45, hjust = 1)
  )

p_all_cov <- ggplot2::ggplot(plot_all, ggplot2::aes(x = bin, y = n_cpg_obs)) +
  ggplot2::geom_area(fill = "#92c5de", alpha = 0.6) +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.5) +
  ggplot2::scale_x_continuous(
    breaks = seq(-meta_window_bp, meta_window_bp, by = 250),
    labels = function(x) paste0(x, " bp")
  ) +
  ggplot2::scale_y_continuous(labels = scales::comma) +
  ggplot2::labs(x = "Distance from silencer center (bp)", y = "CpG observations") +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    axis.text.x      = ggplot2::element_text(angle = 45, hjust = 1)
  )

p_all <- patchwork::wrap_plots(p_all_meta, p_all_cov, ncol = 1, heights = c(3, 1))

ggplot2::ggsave(
  file.path(output_dir,
            sprintf("%s_all_chroms_w%d_b%d_silencer_pbcpg_5mc_metaprofile_combined.pdf",
                    default_sample, meta_window_bp, bin_size)),
  plot = p_all, width = 10, height = 7, device = "pdf"
)
message("Saved combined plot.")

# ── Plot 2: overlay by cCRE class ─────────────────────────────────────────────
message("Plotting overlay by cCRE class...")
plot_class <- meta_class_full[!is.na(fraction_modified_smooth)]

p_overlay_meta <- ggplot2::ggplot(
  plot_class,
  ggplot2::aes(x = bin, y = fraction_modified_smooth, color = ccre_class)
) +
  ggplot2::geom_line(linewidth = 0.7) +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.5) +
  ggplot2::scale_color_manual(values = class_colors, name = "cCRE class") +
  ggplot2::scale_x_continuous(
    breaks = seq(-meta_window_bp, meta_window_bp, by = 250),
    labels = function(x) paste0(x, " bp")
  ) +
  ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  ggplot2::labs(
    title    = "5mCpG Metaprofile at STARR-Silencers by cCRE Class (Overlay)",
    subtitle = sprintf("Sample: %s | All chromosomes | Window: ±%d bp | Bin: %d bp",
                       default_sample, meta_window_bp, bin_size),
    x = "Distance from silencer center (bp)",
    y = "Fraction 5mCpG"
  ) +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::theme(
    plot.title       = ggplot2::element_text(face = "bold", size = 13),
    plot.subtitle    = ggplot2::element_text(size = 10, color = "grey40"),
    panel.grid.minor = ggplot2::element_blank(),
    axis.text.x      = ggplot2::element_text(angle = 45, hjust = 1),
    legend.position  = "right"
  )

p_overlay_cov <- ggplot2::ggplot(
  plot_class,
  ggplot2::aes(x = bin, y = n_cpg_obs, color = ccre_class)
) +
  ggplot2::geom_line(linewidth = 0.5, alpha = 0.7) +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.5) +
  ggplot2::scale_color_manual(values = class_colors, name = "cCRE class") +
  ggplot2::scale_x_continuous(
    breaks = seq(-meta_window_bp, meta_window_bp, by = 250),
    labels = function(x) paste0(x, " bp")
  ) +
  ggplot2::scale_y_continuous(labels = scales::comma) +
  ggplot2::labs(x = "Distance from silencer center (bp)", y = "CpG observations") +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    axis.text.x      = ggplot2::element_text(angle = 45, hjust = 1),
    legend.position  = "right"
  )

p_overlay <- patchwork::wrap_plots(p_overlay_meta, p_overlay_cov, ncol = 1, heights = c(3, 1))

ggplot2::ggsave(
  file.path(output_dir,
            sprintf("%s_all_chroms_w%d_b%d_silencer_pbcpg_5mc_metaprofile_overlay.pdf",
                    default_sample, meta_window_bp, bin_size)),
  plot = p_overlay, width = 11, height = 7, device = "pdf"
)
message("Saved overlay plot.")

# ── Plot 3: faceted by cCRE class ─────────────────────────────────────────────
message("Plotting faceted by cCRE class...")

p_facet_meta <- ggplot2::ggplot(
  plot_class,
  ggplot2::aes(x = bin, y = fraction_modified_smooth, color = ccre_class)
) +
  ggplot2::geom_line(linewidth = 0.7) +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.5) +
  ggplot2::scale_color_manual(values = class_colors, guide = "none") +
  ggplot2::scale_x_continuous(
    breaks = seq(-meta_window_bp, meta_window_bp, by = 500),
    labels = function(x) paste0(x, " bp")
  ) +
  ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  ggplot2::facet_wrap(~ class_label, scales = "free_y") +
  ggplot2::labs(
    title    = "5mCpG Metaprofile at STARR-Silencers by cCRE Class (Faceted)",
    subtitle = sprintf("Sample: %s | All chromosomes | Window: ±%d bp | Bin: %d bp",
                       default_sample, meta_window_bp, bin_size),
    x = "Distance from silencer center (bp)",
    y = "Fraction 5mCpG"
  ) +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(
    plot.title       = ggplot2::element_text(face = "bold", size = 13),
    plot.subtitle    = ggplot2::element_text(size = 10, color = "grey40"),
    panel.grid.minor = ggplot2::element_blank(),
    axis.text.x      = ggplot2::element_text(angle = 45, hjust = 1),
    strip.text       = ggplot2::element_text(face = "bold")
  )

p_facet_cov <- ggplot2::ggplot(
  plot_class,
  ggplot2::aes(x = bin, y = n_cpg_obs, fill = ccre_class)
) +
  ggplot2::geom_area(alpha = 0.5) +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.5) +
  ggplot2::scale_fill_manual(values = class_colors, guide = "none") +
  ggplot2::scale_x_continuous(
    breaks = seq(-meta_window_bp, meta_window_bp, by = 500),
    labels = function(x) paste0(x, " bp")
  ) +
  ggplot2::scale_y_continuous(labels = scales::comma) +
  ggplot2::facet_wrap(~ class_label, scales = "free_y") +
  ggplot2::labs(x = "Distance from silencer center (bp)", y = "CpG obs") +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    axis.text.x      = ggplot2::element_text(angle = 45, hjust = 1),
    strip.text       = ggplot2::element_text(face = "bold")
  )

p_facet <- patchwork::wrap_plots(p_facet_meta, p_facet_cov, ncol = 1, heights = c(3, 1))

ggplot2::ggsave(
  file.path(output_dir,
            sprintf("%s_all_chroms_w%d_b%d_silencer_pbcpg_5mc_metaprofile_facet.pdf",
                    default_sample, meta_window_bp, bin_size)),
  plot = p_facet, width = 14, height = 10, device = "pdf"
)

message("Saved faceted plot.")
message(sprintf("All outputs saved to: %s", output_dir))

 
