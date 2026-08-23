# scripts/tcga/correlation_tcga.R

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

run_tcga_esr1_bub1b <- function(input_csv, output_dir = "results") {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  df <- readr::read_csv(input_csv, show_col_types = FALSE) %>% drop_na(ESR1, BUB1B)
  df <- df %>% filter(sample_type == "Primary Tumor")
  p <- scatter_plot_with_fit(df, "ESR1", "BUB1B", xlab = "ESR1", ylab = "BUB1B", title = "ESR1 vs BUB1B — TCGA")
  ggplot2::ggsave(file.path(output_dir, "tcga_ESR1_BUB1B.png"), p, width = 7, height = 6, dpi = 300)
}
