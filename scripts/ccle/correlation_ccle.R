# scripts/ccle/correlation_ccle.R

# Lightweight wrapper for CCLE correlation analysis (uses R/correlation.R helpers)
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
})

run_ccle_bub1b_mad2l1 <- function(expr_file, meta_file, output_dir = "results") {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  expr <- readr::read_csv(expr_file, show_col_types = FALSE)
  meta <- readr::read_csv(meta_file, show_col_types = FALSE)
  # assume first column is ModelID
  colnames(expr)[1] <- "ModelID"
  bc_ids <- meta %>% filter(OncotreeLineage == "Breast") %>% pull(ModelID)
  df <- expr %>% filter(ModelID %in% bc_ids) %>% select(ModelID, BUB1B, MAD2L1) %>% drop_na() %>% left_join(meta %>% select(ModelID, CellLineName, OncotreeSubtype), by = "ModelID")
  p <- scatter_plot_with_fit(df, "BUB1B", "MAD2L1", xlab = "BUB1B", ylab = "MAD2L1", title = "BUB1B vs MAD2L1 — CCLE")
  ggplot2::ggsave(file.path(output_dir, "ccle_bub1b_mad2l1.png"), p, width = 7, height = 6, dpi = 300)
  ggplot2::ggsave(file.path(output_dir, "ccle_bub1b_mad2l1.pdf"), p, width = 7, height = 6)
  write.csv(df %>% select(CellLineName, OncotreeSubtype, BUB1B, MAD2L1), file.path(output_dir, "ccle_celllines_used.csv"), row.names = FALSE)
}
