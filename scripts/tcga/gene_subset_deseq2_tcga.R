# scripts/tcga/gene_subset_deseq2_tcga.R

# Wrapper script demonstrating use of DESeq2 on a gene subset for a TCGA SummarizedExperiment
# This is a dataset-specific example that calls R/expression.R helpers.

suppressPackageStartupMessages({
  library(DESeq2)
  library(readr)
  library(SummarizedExperiment)
})

run_gene_subset_deseq2 <- function(counts_file, genes_file, output_dir = "results") {
  stopifnot(file.exists(counts_file))
  stopifnot(file.exists(genes_file))
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  counts_data <- readRDS(counts_file)
  genes_list <- readr::read_csv(genes_file, show_col_types = FALSE)
  stopifnot(inherits(counts_data, "SummarizedExperiment"))
  count_matrix <- SummarizedExperiment::assay(counts_data)
  sample_info <- as.data.frame(SummarizedExperiment::colData(counts_data))
  rownames(count_matrix) <- sub("\\..*$", "", rownames(count_matrix))
  ## Infer condition from TCGA barcode if present
  if ("barcode" %in% colnames(sample_info)) {
    sample_code <- sapply(strsplit(as.character(sample_info$barcode), "-"), `[`, 4)
    sample_info$condition <- ifelse(grepl("^01", sample_code), "Tumor", ifelse(grepl("^11", sample_code), "Normal", NA))
  }
  count_matrix <- count_matrix[, rownames(sample_info), drop = FALSE]
  genes_of_interest <- genes_list$ENSEMBL
  count_matrix <- count_matrix[rownames(count_matrix) %in% genes_of_interest, , drop = FALSE]
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = round(count_matrix), colData = sample_info, design = ~ condition)
  dds <- DESeq2::DESeq(dds)
  res <- DESeq2::results(dds)
  # write results
  out_file <- file.path(output_dir, paste0("gene_subset_DESeq2_results.csv"))
  write.csv(data.frame(ENSEMBL = rownames(res), LOG2FC = res$log2FoldChange, p_value = res$pvalue, adjusted_p_value = res$padj), out_file, row.names = FALSE)
  message("Wrote: ", out_file)
}
