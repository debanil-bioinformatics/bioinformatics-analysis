############################################################
# Author: Debanil Dhar
# Year: 2025
#
# Description:
# Aggregates gene-level differential expression results
# across multiple TCGA projects by computing a composite
# score based on effect size (log2FC) and statistical
# significance (adjusted p-value), followed by rank
# aggregation across datasets.
#
# Input:
# - Multiple CSV files containing differential expression
#   results with log2 fold change and adjusted p-values
#
# Output:
# - A ranked summary table of genes across all datasets
#
# Notes:
# - Designed to integrate results from multiple TCGA
#   cancer types or analysis batches
############################################################

## Load required libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
})

## ---------------- USER PARAMETERS ----------------
input_dir   <- "lfc_results"       # Directory with DE result CSVs
file_regex  <- "^lfc_results_.*\\.csv$"
output_file <- "tcga_gene_rank_summary.csv"
## -------------------------------------------------

## Sanity check
stopifnot(dir.exists(input_dir))

## List input files
csv_files <- list.files(
  path    = input_dir,
  pattern = file_regex,
  full.names = TRUE
)

stopifnot(length(csv_files) > 0)

## Scoring function
calculate_score <- function(log2fc, padj) {
  padj <- pmax(padj, 1e-300)  # Avoid -log10(0)
  abs(log2fc) * -log10(padj)
}

## Process a single result file
process_file <- function(file) {
  
  df <- read.csv(file, stringsAsFactors = FALSE)
  
  required_cols <- c("ENSEMBL", "SYMBOL", "log2FoldChange", "padj")
  if (!all(required_cols %in% colnames(df))) {
    stop("Missing required columns in: ", basename(file))
  }
  
  df %>%
    filter(!is.na(log2FoldChange), !is.na(padj)) %>%
    mutate(
      score = calculate_score(log2FoldChange, padj),
      dataset = gsub("^lfc_results_|\\.csv$", "", basename(file))
    ) %>%
    arrange(desc(score)) %>%
    mutate(rank = row_number()) %>%
    select(ENSEMBL, SYMBOL, score, rank, dataset)
}

## Process all files
all_data <- map_dfr(csv_files, process_file)

## Aggregate rankings across datasets
rank_summary <- all_data %>%
  group_by(ENSEMBL, SYMBOL) %>%
  summarise(
    mean_score      = mean(score, na.rm = TRUE),
    median_rank     = median(rank, na.rm = TRUE),
    mean_rank       = mean(rank, na.rm = TRUE),
    datasets_present = n(),
    .groups = "drop"
  ) %>%
  arrange(median_rank, mean_rank) %>%
  mutate(overall_rank = row_number())

## Write output
write.csv(rank_summary, output_file, row.names = FALSE)

message("Rank aggregation completed successfully.")
message("Top 10 ranked genes:")
print(head(rank_summary, 10))
