source("/project/spott/cshan/fiber-seq/plot_nuc_distribution_functions.R")


cage_tss_gene_path <- "/project/spott/cshan/annotations/fantom5/fantom5.hg38.LCL.consensus.CAGE_peaks.withGene.bed.gz"
sample_root_all_chr <- "/project/spott/1_Shared_projects/LCL_Fiber_seq/FIRE/results"
default_sample_all_chr <- "AL10_bc2178_19130"
window_bp_all_chr <- 2000L

read_nuc_bed12 <- function(path) {
  fread(path, header = FALSE,
        col.names = c("chrom","chromStart","chromEnd","name","score","strand",
                      "thickStart","thickEnd","itemRgb",
                      "blockCount","blockSizes","blockStarts"))
}

load_all_chr_nucleosome_midpoints <- function(bed_files) {
  nucs_all <- rbindlist(lapply(bed_files, function(path) {
    calculate_nucleosome_midpoints(read_nuc_bed12(path))
  }), use.names = TRUE, fill = TRUE)
  nucs_all[order(chrom, nuc_mid)]
}

granges_to_dt <- function(gr) {
  data.table(
    chrom = as.character(seqnames(gr)),
    start = start(gr) - 1L,
    end = end(gr),
    strand = as.character(strand(gr)),
    gene_id = as.character(mcols(gr)$gene_id),
    score = as.character(mcols(gr)$score)
  )
}

find_nucs_near_tss_with_gene <- function(nuc_gr, tss_gr, window_bp = 2000L) {
  tss_windows <- promoters(tss_gr, upstream = window_bp, downstream = window_bp + 1L)
  hits <- findOverlaps(nuc_gr, tss_windows)
  if (!length(hits)) {
    return(data.table(
      gene_id = character(),
      gene_id_base = character(),
      seqnames = character(),
      strand = character(),
      tss_index = integer(),
      nuc_pos = integer(),
      tss_pos = integer(),
      rel_pos = integer()
    ))
  }

  query_idx <- queryHits(hits)
  subject_idx <- subjectHits(hits)
  nuc_pos <- start(nuc_gr)[query_idx]
  tss_pos <- start(tss_gr)[subject_idx]
  tss_str <- as.character(strand(tss_gr))[subject_idx]
  gene_id <- as.character(mcols(tss_gr)$gene_id[subject_idx])

  data.table(
    gene_id = gene_id,
    gene_id_base = sub("\\..*$", "", gene_id),
    seqnames = as.character(seqnames(tss_gr))[subject_idx],
    strand = tss_str,
    tss_index = subject_idx,
    nuc_pos = nuc_pos,
    tss_pos = tss_pos,
    rel_pos = ifelse(tss_str == "+", nuc_pos - tss_pos, tss_pos - nuc_pos)
  )
}

summarize_gene_tss_distance <- function(hit_dt) {
  if (!nrow(hit_dt)) {
    return(data.table(
      gene_id = character(),
      gene_id_base = character(),
      mean_rel_pos_bp = numeric(),
      median_rel_pos_bp = numeric(),
      mean_abs_rel_pos_bp = numeric(),
      median_abs_rel_pos_bp = numeric(),
      sd_rel_pos_bp = numeric(),
      n_midpoints = integer(),
      n_tss = integer(),
      n_chromosomes = integer()
    ))
  }

  hit_dt[, .(
    mean_rel_pos_bp = as.numeric(mean(rel_pos)),
    median_rel_pos_bp = as.numeric(median(rel_pos)),
    mean_abs_rel_pos_bp = as.numeric(mean(abs(rel_pos))),
    median_abs_rel_pos_bp = as.numeric(median(abs(rel_pos))),
    sd_rel_pos_bp = as.numeric(stats::sd(rel_pos)),
    n_midpoints = .N,
    n_tss = uniqueN(tss_index),
    n_chromosomes = uniqueN(seqnames)
  ), by = .(gene_id, gene_id_base)][order(mean_abs_rel_pos_bp, -n_midpoints, gene_id)]
}

selected_sample_all_chr <- file.path(sample_root_all_chr, default_sample_all_chr)
nuc_by_chr_dir_all_chr <- file.path(selected_sample_all_chr, "extracted_results/nuc_by_chr")
bed_files_all_chr <- sort(list.files(
  nuc_by_chr_dir_all_chr,
  pattern = "ft_extracted_nuc\\.chr.*\\.bed\\.gz$",
  full.names = TRUE
))
if (!length(bed_files_all_chr)) {
  stop("No chromosome-level nucleosome BED files found in ", nuc_by_chr_dir_all_chr)
}

tss_gr_all_chr <- load_tss_annotations_bed9_with_gene_id(cage_tss_gene_path)
tss_gr_all_chr <- tss_gr_all_chr[
  !is.na(mcols(tss_gr_all_chr)$gene_id) &
    mcols(tss_gr_all_chr)$gene_id != "NA"
]

nucs_all_chr <- load_all_chr_nucleosome_midpoints(bed_files_all_chr)
nuc_gr_all_chr <- nuc_midpoints_to_granges(nucs_all_chr)
cage_gene_hits_all_chr <- find_nucs_near_tss_with_gene(
  nuc_gr_all_chr,
  tss_gr_all_chr,
  window_bp = window_bp_all_chr
)
if (!nrow(cage_gene_hits_all_chr)) {
  stop("No nucleosome midpoints overlapped the CAGE TSS windows across all chromosomes.")
}

meta_raw_all_chr <- make_position_meta(
  nuc_gr_all_chr,
  tss_gr_all_chr,
  window_bp = window_bp_all_chr,
  bin_size = 10L,
  normalize_by_tss = FALSE,
  label = "CAGE TSS"
)
meta_norm_all_chr <- make_position_meta(
  nuc_gr_all_chr,
  tss_gr_all_chr,
  window_bp = window_bp_all_chr,
  bin_size = 10L,
  normalize_by_tss = TRUE,
  label = "CAGE TSS"
)

cage_gene_distance_summary <- summarize_gene_tss_distance(cage_gene_hits_all_chr)
print(head(cage_gene_distance_summary, 20))

meta_plot_all_chr <- ggplot(meta_raw_all_chr, aes(x = bin, y = N)) +
  geom_line(linewidth = 0.6) +
  geom_smooth(method = "loess", span = 0.05, se = FALSE,
              color = "steelblue", linewidth = 1) +
  geom_vline(xintercept = 0, linetype = "dashed",
             color = "red", linewidth = 0.5) +
  scale_x_continuous(breaks = seq(-window_bp_all_chr, window_bp_all_chr, 500),
                     labels = function(x) paste0(x / 1000, " kb")) +
  labs(
    x = "Distance to TSS",
    y = "Nucleosome midpoint count",
    title = "Nucleosome midpoint distribution around CAGE TSS",
    subtitle = sprintf(
      "Sample: %s | all chromosomes | n nucleosomes = %s | n TSS = %s",
      default_sample_all_chr,
      format(length(nuc_gr_all_chr), big.mark = ","),
      format(length(tss_gr_all_chr), big.mark = ",")
    )
  ) +
  theme_bw(base_size = 12)
print(meta_plot_all_chr)

meta_plot_norm_all_chr <- ggplot(meta_norm_all_chr, aes(x = bin, y = N)) +
  geom_line(linewidth = 0.7, color = "#d62728") +
  geom_smooth(method = "loess", span = 0.06, se = FALSE,
              linewidth = 1, color = "#d62728") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red", linewidth = 0.5) +
  scale_x_continuous(breaks = seq(-window_bp_all_chr, window_bp_all_chr, 500),
                     labels = function(x) paste0(x / 1000, " kb")) +
  labs(
    x = "Distance to TSS",
    y = "Avg nucleosome midpoints per TSS (10 bp bins)",
    title = "CAGE TSS",
    subtitle = sprintf("Sample: %s | all chromosomes", default_sample_all_chr)
  ) +
  theme_bw(base_size = 12)
print(meta_plot_norm_all_chr)

gene_distance_plot <- ggplot(cage_gene_distance_summary, aes(x = mean_abs_rel_pos_bp)) +
  geom_histogram(binwidth = 25, fill = "#be1c54ff", color = "white") +
  labs(
    x = "Mean absolute distance to TSS per gene (bp)",
    y = "Number of genes",
    title = "Per-gene average nucleosome distance to CAGE TSS",
    subtitle = sprintf(
      "Sample: %s | all chromosomes | genes = %s",
      default_sample_all_chr,
      format(nrow(cage_gene_distance_summary), big.mark = ",")
    )
  ) +
  theme_bw(base_size = 12)
print(gene_distance_plot)

out_dir <- "/project/spott/cshan/fiber-seq/results/plots/nuc"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

fwrite(
  data.table(
    sample = default_sample_all_chr,
    chrom_bed = basename(bed_files_all_chr),
    chrom_bed_path = bed_files_all_chr
  ),
  file.path(out_dir,
            sprintf("%s_CAGE_TSS_input_bed_files_all_chr.tsv",
                    default_sample_all_chr)),
  sep = "\t"
)

fwrite(
  nucs_all_chr,
  file.path(out_dir,
            sprintf("%s_CAGE_TSS_nucleosome_midpoints_all_chr.tsv",
                    default_sample_all_chr)),
  sep = "\t"
)

fwrite(
  granges_to_dt(tss_gr_all_chr),
  file.path(out_dir,
            sprintf("%s_CAGE_TSS_annotations_all_chr.tsv",
                    default_sample_all_chr)),
  sep = "\t"
)

fwrite(
  cage_gene_hits_all_chr,
  file.path(out_dir,
            sprintf("%s_CAGE_TSS_nucleosome_tss_hits_all_chr.tsv",
                    default_sample_all_chr)),
  sep = "\t"
)

fwrite(
  meta_raw_all_chr,
  file.path(out_dir,
            sprintf("%s_CAGE_TSS_nucleosome_midpoints_raw_all_chr.tsv",
                    default_sample_all_chr)),
  sep = "\t"
)

fwrite(
  meta_norm_all_chr,
  file.path(out_dir,
            sprintf("%s_CAGE_TSS_nucleosome_midpoints_perTSS_all_chr.tsv",
                    default_sample_all_chr)),
  sep = "\t"
)

fwrite(
  cage_gene_distance_summary,
  file.path(out_dir,
            sprintf("%s_CAGE_TSS_gene_mean_distance_all_chr.tsv",
                    default_sample_all_chr)),
  sep = "\t"
)

ggsave(
  file.path(out_dir,
            sprintf("%s_CAGE_TSS_nucleosome_midpoints_raw_all_chr.pdf",
                    default_sample_all_chr)),
  plot = meta_plot_all_chr,
  width = 6,
  height = 4
)

ggsave(
  file.path(out_dir,
            sprintf("%s_CAGE_TSS_nucleosome_midpoints_perTSS_all_chr.pdf",
                    default_sample_all_chr)),
  plot = meta_plot_norm_all_chr,
  width = 6,
  height = 4
)

ggsave(
  file.path(out_dir,
            sprintf("%s_CAGE_TSS_gene_mean_abs_distance_all_chr.pdf",
                    default_sample_all_chr)),
  plot = gene_distance_plot,
  width = 7,
  height = 4.5
)
