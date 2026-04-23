source("/project/spott/cshan/fiber-seq/code/footprintR_modbam_code/modbam_footprintR_functions.R")

parse_cli_args <- function(argv) {
  out <- list(
    sample = NULL,
    chrom = NULL,
    start = NULL,
    end = NULL,
    region = NULL,
    metadata = "/project/spott/1_Shared_projects/LCL_Fiber_seq/Data/LCL_sample_metatable_merged_samples_31samples.csv",
    output_root = "/project/spott/cshan/fiber-seq/results/region_m6a_footprint_overlay",
    chrom_sizes = "/project/spott/cshan/annotations/hg38.chrom.sizes",
    mod_prob_threshold = 0.9
  )

  i <- 1L
  while (i <= length(argv)) {
    arg <- argv[[i]]
    if (grepl("^--[^=]+=", arg)) {
      pieces <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1]]
      key <- pieces[[1]]
      value <- paste(pieces[-1], collapse = "=")
    } else if (grepl("^--", arg)) {
      key <- sub("^--", "", arg)
      if (i == length(argv)) {
        stop(sprintf("Missing value for argument: %s", arg), call. = FALSE)
      }
      i <- i + 1L
      value <- argv[[i]]
    } else if (grepl("=", arg, fixed = TRUE)) {
      pieces <- strsplit(arg, "=", fixed = TRUE)[[1]]
      key <- pieces[[1]]
      value <- paste(pieces[-1], collapse = "=")
    } else {
      stop(sprintf("Unrecognized argument: %s", arg), call. = FALSE)
    }

    out[[key]] <- value
    i <- i + 1L
  }

  out
}

parse_region_string <- function(region) {
  m <- regexec("^([^:]+):(\\d+)-(\\d+)$", region)
  parts <- regmatches(region, m)[[1]]
  if (length(parts) != 4L) {
    stop("Region must look like chr1:10000-20000", call. = FALSE)
  }

  list(
    chrom = parts[[2]],
    start = as.integer(parts[[3]]),
    end = as.integer(parts[[4]])
  )
}

normalize_chrom <- function(chrom) {
  chrom <- as.character(chrom)
  if (startsWith(chrom, "chr")) chrom else paste0("chr", chrom)
}

args <- parse_cli_args(commandArgs(trailingOnly = TRUE))

if (!is.null(args$region) && (is.null(args$chrom) || is.null(args$start) || is.null(args$end))) {
  region_parts <- parse_region_string(args$region)
  args$chrom <- region_parts$chrom
  args$start <- region_parts$start
  args$end <- region_parts$end
}

if (is.null(args$sample) || is.null(args$chrom) || is.null(args$start) || is.null(args$end)) {
  stop("Provide --sample plus either --region or --chrom/--start/--end.", call. = FALSE)
}

sample_name       <- as.character(args$sample)
chrom             <- normalize_chrom(args$chrom)
region_start      <- as.integer(args$start)
region_end        <- as.integer(args$end)
mod_prob_threshold <- as.numeric(args$mod_prob_threshold)
chrom_sizes_path  <- as.character(args$chrom_sizes)

if (is.na(region_start) || is.na(region_end) || region_start < 1L || region_end < region_start) {
  stop("Invalid region coordinates.", call. = FALSE)
}

if (!file.exists(chrom_sizes_path)) {
  stop(sprintf("chrom sizes file not found: %s", chrom_sizes_path), call. = FALSE)
}

metadata_path <- args$metadata
output_root   <- args$output_root
region_slug   <- sprintf("%s_%d_%d", chrom, region_start, region_end)
region_dir    <- file.path(output_root, sample_name, chrom, region_slug)
dir.create(region_dir, recursive = TRUE, showWarnings = FALSE)

meta <- data.table::fread(metadata_path)
sample_row <- meta[meta$sample_name == sample_name, ]
if (!nrow(sample_row)) {
  stop(sprintf("Sample %s was not found in %s", sample_name, metadata_path), call. = FALSE)
}
sample_row <- sample_row[1]

bam_path <- sample_row$bam_file[[1]]
assert_file_exists(bam_path, "sample BAM")
assert_file_exists(paste0(bam_path, ".bai"), "sample BAM index")

bam_seqinfo <- read_bam_seqinfo(bam_path)
if (!(chrom %in% GenomeInfoDb::seqnames(bam_seqinfo))) {
  stop(sprintf("Chromosome %s is absent from %s", chrom, bam_path), call. = FALSE)
}

region_gr <- GenomicRanges::GRanges(
  seqnames = chrom,
  ranges = IRanges::IRanges(start = region_start, end = region_end)
)
GenomeInfoDb::seqinfo(region_gr) <- bam_seqinfo[chrom]
region_gr <- GenomicRanges::trim(region_gr)

region_info <- data.table::data.table(
  sample_name = sample_name,
  chromosome  = chrom,
  region_start = start(region_gr),
  region_end   = end(region_gr),
  bam_path     = bam_path,
  metadata_path = metadata_path
)
write_gz_tsv(region_info, file.path(region_dir, sprintf("%s.%s.region_info.tsv.gz", sample_name, region_slug)))

position_dt <- read_modbam_summary_positions(
  bamfiles          = setNames(bam_path, sample_name),
  sample_name       = sample_name,
  promoter_windows  = region_gr,
  modbase_code      = "a",
  seqinfo           = bam_seqinfo,
  bpparam           = BiocParallel::SerialParam(),
  modProbThreshold  = mod_prob_threshold
)

position_tsv <- file.path(region_dir, sprintf("%s.%s.m6A_position_calls.tsv.gz", sample_name, region_slug))
write_gz_tsv(position_dt, position_tsv)

bedgraph_path <- file.path(region_dir, sprintf("%s.%s.m6A_fraction.bedgraph", sample_name, region_slug))
if (nrow(position_dt)) {
  bedgraph_dt <- position_dt[, .(
    score = mean(frac_modified, na.rm = TRUE)
  ), by = .(chrom = call_chromosome, start0 = position - 1L, end = position)]
  bedgraph_dt <- bedgraph_dt[order(chrom, start0, end)]
} else {
  bedgraph_dt <- data.table::data.table(
    chrom  = character(),
    start0 = integer(),
    end    = integer(),
    score  = numeric()
  )
}
data.table::fwrite(bedgraph_dt, file = bedgraph_path, sep = "\t", col.names = FALSE)

# Convert bedGraph to bigWig using the genome-wide hg38.chrom.sizes.
# The old approach wrote a per-run chrom sizes file derived from the BAM
# seqinfo, which often had missing or incorrect chromosome lengths, causing
# bedGraphToBigWig to exit with status 255.
bigwig_path <- file.path(region_dir, sprintf("%s.%s.m6A_fraction.bw", sample_name, region_slug))

bedgraph_to_bigwig <- Sys.which("bedGraphToBigWig")
if (!nzchar(bedgraph_to_bigwig) && file.exists("/project/spott/cshan/tools/bedGraphToBigWig")) {
  bedgraph_to_bigwig <- "/project/spott/cshan/tools/bedGraphToBigWig"
}

if (nzchar(bedgraph_to_bigwig) && nrow(bedgraph_dt)) {
  conv_result <- system2(
    command = bedgraph_to_bigwig,
    args    = c(bedgraph_path, chrom_sizes_path, bigwig_path),
    stdout  = TRUE,
    stderr  = TRUE
  )
  conv_status <- attr(conv_result, "status")
  if (!is.null(conv_status) && conv_status != 0L) {
    warning(sprintf(
      "bedGraphToBigWig exited with status %d for %s. bigWig will not be available.\nOutput:\n%s",
      conv_status, bigwig_path, paste(conv_result, collapse = "\n")
    ))
    # Remove the broken output file so downstream code does not try to open it
    if (file.exists(bigwig_path)) file.remove(bigwig_path)
    bigwig_path <- NULL
  }
} else {
  bigwig_path <- NULL
}

plot_path <- file.path(region_dir, sprintf("%s.%s.m6A_fraction.pdf", sample_name, region_slug))
if (nrow(position_dt)) {
  p <- ggplot2::ggplot(position_dt, ggplot2::aes(x = position, y = frac_modified)) +
    ggplot2::geom_line(color = "#7B3294", linewidth = 0.5) +
    ggplot2::coord_cartesian(xlim = c(start(region_gr), end(region_gr)), ylim = c(0, 1)) +
    ggplot2::labs(
      title = sprintf("Region m6A fraction: %s %s", sample_name, region_slug),
      x     = sprintf("Genomic position on %s", chrom),
      y     = "Fraction modified"
    ) +
    ggplot2::theme_bw(base_size = 11)
} else {
  p <- ggplot2::ggplot() +
    ggplot2::annotate("text", x = 0, y = 0, label = "No m6A positions with coverage in this region.") +
    ggplot2::theme_void() +
    ggplot2::labs(title = sprintf("Region m6A fraction: %s %s", sample_name, region_slug))
}
ggplot2::ggsave(plot_path, p, width = 8, height = 2.8)

summary_path <- file.path(region_dir, sprintf("%s.%s.m6A_region_summary.tsv.gz", sample_name, region_slug))
summary_dt <- if (nrow(position_dt)) {
  position_dt[, .(
    modified_calls          = sum(modified_calls, na.rm = TRUE),
    total_calls             = sum(total_calls, na.rm = TRUE),
    mean_fraction_modified  = mean(frac_modified, na.rm = TRUE),
    median_fraction_modified = stats::median(frac_modified, na.rm = TRUE),
    covered_positions       = .N
  )]
} else {
  data.table::data.table(
    modified_calls          = 0L,
    total_calls             = 0L,
    mean_fraction_modified  = NA_real_,
    median_fraction_modified = NA_real_,
    covered_positions       = 0L
  )
}
summary_dt[, `:=`(
  sample_name        = sample_name,
  chromosome         = chrom,
  region_start       = start(region_gr),
  region_end         = end(region_gr),
  mod_prob_threshold = mod_prob_threshold
)]
data.table::setcolorder(
  summary_dt,
  c(
    "sample_name", "chromosome", "region_start", "region_end", "mod_prob_threshold",
    "modified_calls", "total_calls", "mean_fraction_modified", "median_fraction_modified",
    "covered_positions"
  )
)
write_gz_tsv(summary_dt, summary_path)

message("Region m6A outputs written to: ", region_dir)

# ── 5mC (methylated-C) section ────────────────────────────────────────────────
position_dt_5mc <- read_modbam_summary_positions(
  bamfiles          = setNames(bam_path, sample_name),
  sample_name       = sample_name,
  promoter_windows  = region_gr,
  modbase_code      = "m",
  seqinfo           = bam_seqinfo,
  bpparam           = BiocParallel::SerialParam(),
  modProbThreshold  = mod_prob_threshold
)

bedgraph_path_5mc <- file.path(region_dir, sprintf("%s.%s.5mC_fraction.bedgraph", sample_name, region_slug))
if (nrow(position_dt_5mc)) {
  bedgraph_dt_5mc <- position_dt_5mc[, .(
    score = mean(frac_modified, na.rm = TRUE)
  ), by = .(chrom = call_chromosome, start0 = position - 1L, end = position)]
  bedgraph_dt_5mc <- bedgraph_dt_5mc[order(chrom, start0, end)]
} else {
  bedgraph_dt_5mc <- data.table::data.table(
    chrom  = character(),
    start0 = integer(),
    end    = integer(),
    score  = numeric()
  )
}
data.table::fwrite(bedgraph_dt_5mc, file = bedgraph_path_5mc, sep = "\t", col.names = FALSE)

bigwig_path_5mc <- file.path(region_dir, sprintf("%s.%s.5mC_fraction.bw", sample_name, region_slug))
if (nzchar(bedgraph_to_bigwig) && nrow(bedgraph_dt_5mc)) {
  conv_result_5mc <- system2(
    command = bedgraph_to_bigwig,
    args    = c(bedgraph_path_5mc, chrom_sizes_path, bigwig_path_5mc),
    stdout  = TRUE,
    stderr  = TRUE
  )
  conv_status_5mc <- attr(conv_result_5mc, "status")
  if (!is.null(conv_status_5mc) && conv_status_5mc != 0L) {
    warning(sprintf(
      "bedGraphToBigWig exited with status %d for %s.\nOutput:\n%s",
      conv_status_5mc, bigwig_path_5mc, paste(conv_result_5mc, collapse = "\n")
    ))
    if (file.exists(bigwig_path_5mc)) file.remove(bigwig_path_5mc)
    bigwig_path_5mc <- NULL
  }
} else {
  bigwig_path_5mc <- NULL
}

# Print the bigWig paths (or empty strings) so the calling pipeline can pick them up
cat(if (!is.null(bigwig_path)) bigwig_path else "", "\n", sep = "")
cat(if (!is.null(bigwig_path_5mc)) bigwig_path_5mc else "", "\n", sep = "")