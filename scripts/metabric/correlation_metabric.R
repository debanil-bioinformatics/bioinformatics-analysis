# scripts/metabric/correlation_metabric.R

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

run_metabric_correlation <- function(input_csv, output_dir = "results") {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  df <- readr::read_csv(input_csv, show_col_types = FALSE)
  p <- scatter_plot_with_fit(df, "GeneX", "GeneY", title = "METABRIC correlation example")
  ggplot2::ggsave(file.path(output_dir, "metabric_corr.png"), p, width = 7, height = 6)
}
