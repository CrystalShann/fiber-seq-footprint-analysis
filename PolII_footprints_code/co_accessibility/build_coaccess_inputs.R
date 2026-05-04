#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(GenomicRanges))

args <- commandArgs(trailingOnly = TRUE)
n_genes <- if (length(args) >= 1) as.integer(args[[1]]) else 500L
if (is.na(n_genes) || n_genes <= 0) {
  stop("First argument must be a positive integer number of high/low genes to write.")
}

pausing_index_file <- "/project/spott/cshan/fiber-seq/results/PolII/annotations/pausing_index_principal_with_CAGE_TSS.csv"
cre_bed <- "/project/spott/cshan/annotations/GRCh38-cCREs.bed"
out_dir <- "/project/spott/cshan/fiber-seq/results/co_accessibility/inputs"

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

principal_with_cage <- data.table::fread(pausing_index_file)
required_cols <- c("gene_name", "gene_id", "chrom", "strand", "gstart", "gend", "tss", "PI")
missing_cols <- setdiff(required_cols, colnames(principal_with_cage))
if (length(missing_cols) > 0) {
  stop("Missing required pausing-index columns: ", paste(missing_cols, collapse = ", "))
}

genes_df <- principal_with_cage %>%
  dplyr::mutate(
    PI = as.numeric(PI),
    cage_tss = as.numeric(cage_tss),
    tss = as.numeric(tss),
    tss_for_window = dplyr::if_else(!is.na(cage_tss), cage_tss, tss),
    start = pmin(as.integer(gstart), as.integer(gend)),
    end = pmax(as.integer(gstart), as.integer(gend))
  ) %>%
  dplyr::filter(
    !is.na(gene_name),
    gene_name != "",
    !is.na(chrom),
    !is.na(start),
    !is.na(end),
    !is.na(tss_for_window),
    is.finite(PI)
  ) %>%
  dplyr::arrange(gene_name, dplyr::desc(!is.na(cage_tss)), dplyr::desc(PI)) %>%
  dplyr::distinct(gene_name, .keep_all = TRUE) %>%
  dplyr::mutate(
    annot_start = pmin(start, as.integer(tss_for_window)),
    annot_end = pmax(end, as.integer(tss_for_window))
  )

genes_gr <- GenomicRanges::makeGRangesFromDataFrame(
  genes_df %>%
    dplyr::transmute(
      seqnames = chrom,
      start = annot_start,
      end = annot_end,
      strand = strand,
      gene_name = gene_name,
      gene_id = gene_id,
      transcript_id = transcript_id,
      PI = PI,
      tss = as.integer(tss_for_window),
      original_tss = as.integer(tss),
      cage_tss = as.integer(cage_tss),
      gstart = as.integer(gstart),
      gend = as.integer(gend)
    ),
  keep.extra.columns = TRUE
)

cre_df <- data.table::fread(cre_bed)
if (ncol(cre_df) < 6) {
  stop("Expected at least 6 columns in cCRE BED: ", cre_bed)
}
cre_df <- cre_df[, 1:6]
colnames(cre_df) <- c("chr", "start", "end", "accession1", "accession2", "CRE_label")
cre_df <- cre_df %>%
  dplyr::mutate(
    start = as.integer(start) + 1L,
    end = as.integer(end),
    CRE_ID = paste0(accession1, ".", accession2, ".", CRE_label)
  )
CREs_gr <- GenomicRanges::makeGRangesFromDataFrame(cre_df, keep.extra.columns = TRUE)

genomic_annots <- list(
  genes = genes_gr,
  CREs = CREs_gr,
  metadata = list(
    genome = "hg38",
    gene_source = pausing_index_file,
    cre_source = cre_bed,
    tss_preference = "cage_tss, falling back to tss",
    created = as.character(Sys.time())
  )
)

genomic_annots_rds <- file.path(out_dir, "genomic.annots.hg38.rds")
saveRDS(genomic_annots, genomic_annots_rds)

ranked_genes <- genes_df %>%
  dplyr::arrange(dplyr::desc(PI), gene_name)

high_genes <- ranked_genes %>%
  dplyr::slice_head(n = n_genes) %>%
  dplyr::pull(gene_name)

low_genes <- ranked_genes %>%
  dplyr::arrange(PI, gene_name) %>%
  dplyr::slice_head(n = n_genes) %>%
  dplyr::pull(gene_name)

high_gene_list_file <- file.path(out_dir, "high_pausing_genes.txt")
low_gene_list_file <- file.path(out_dir, "low_pausing_genes.txt")
data.table::fwrite(data.table::data.table(gene_name = high_genes), high_gene_list_file, col.names = FALSE)
data.table::fwrite(data.table::data.table(gene_name = low_genes), low_gene_list_file, col.names = FALSE)

ranked_file <- file.path(out_dir, "pausing_index_genes_ranked_with_CAGE_TSS.tsv")
data.table::fwrite(
  genes_df %>%
    dplyr::arrange(dplyr::desc(PI), gene_name) %>%
    dplyr::select(gene_name, gene_id, transcript_id, chrom, strand, gstart, gend, tss, cage_tss, tss_for_window, PI),
  ranked_file,
  sep = "\t"
)

message("Wrote: ", genomic_annots_rds)
message("Wrote: ", high_gene_list_file)
message("Wrote: ", low_gene_list_file)
message("Wrote: ", ranked_file)
