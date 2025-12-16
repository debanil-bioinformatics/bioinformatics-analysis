############################################################
# Author: Debanil Dhar
# Year: 2025
#
# Description:
# Generalized differential expression analysis of a
# user-defined gene subset between tumor and normal
# samples across TCGA projects using DESeq2.
#
# Input:
# - TCGA STAR-count SummarizedExperiment (RDS)
# - CSV file containing ENSEMBL gene IDs of interest
#
# Output:
# - Project-specific CSV files with log2 fold change,
#   p-values, and adjusted p-values
#
# Notes:
# - Tumor (01) vs Normal (11) samples are inferred from
#   TCGA barcodes
# - Designed to be reusable across TCGA cancer types
############################################################

## Load required libraries
suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(readr)
  library(SummarizedExperiment)
  library(biomaRt)
})

## ---------------- USER PARAMETERS ----------------
tcga_project <- "TCGA-ESCA"      # e.g., TCGA-BRCA, TCGA-LUAD
genes_file   <- "gene_list.csv" # CSV with ENSEMBL column
## -------------------------------------------------

## Define input/output paths
counts_file <- paste0(tcga_project, "_STAR_counts.rds")
output_dir  <- file.path("results", tcga_project)

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

## Sanity checks
stopifnot(file.exists(counts_file))
stopifnot(file.exists(genes_file))

## Load data
counts_data <- readRDS(counts_file)
genes_list  <- read_csv(genes_file, show_col_types = FALSE)

stopifnot(inherits(counts_data, "SummarizedExperiment"))

## Extract count matrix and sample metadata
count_matrix <- assay(counts_data)
sample_info  <- as.data.frame(colData(counts_data))

## Remove Ensembl version numbers
rownames(count_matrix) <- sub("\\..*$", "", rownames(count_matrix))

## Parse TCGA barcodes
sample_info <- sample_info %>%
  mutate(
    sample_code = sapply(strsplit(as.character(barcode), "-"), `[`, 4),
    condition = case_when(
      grepl("^01", sample_code) ~ "Tumor",
      grepl("^11", sample_code) ~ "Normal",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(condition))

## Subset count matrix
count_matrix <- count_matrix[, rownames(sample_info), drop = FALSE]

## Subset to genes of interest
genes_of_interest <- genes_list$ENSEMBL
count_matrix <- count_matrix[rownames(count_matrix) %in% genes_of_interest, , drop = FALSE]

## Map Ensembl IDs to gene symbols
mart <- useMart(
  biomart = "ensembl",
  dataset = "hsapiens_gene_ensembl",
  host    = "https://useast.ensembl.org"
)

gene_map <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters    = "ensembl_gene_id",
  values     = rownames(count_matrix),
  mart       = mart
)

## Differential expression analysis
dds <- DESeqDataSetFromMatrix(
  countData = round(count_matrix),
  colData   = sample_info,
  design    = ~ condition
)

dds <- DESeq(dds)
res <- results(dds)

## Assemble results
results_df <- data.frame(
  ENSEMBL = rownames(res),
  SYMBOL  = gene_map$hgnc_symbol[match(rownames(res), gene_map$ensembl_gene_id)],
  LOG2FC  = res$log2FoldChange,
  p_value = res$pvalue,
  adjusted_p_value = res$padj
) %>%
  filter(!is.na(SYMBOL)) %>%
  arrange(desc(LOG2FC))

## Write output
output_file <- file.path(
  output_dir,
  paste0(tcga_project, "_gene_subset_DESeq2_results.csv")
)

write.csv(results_df, output_file, row.names = FALSE)

message("DESeq2 analysis completed for ", tcga_project)
