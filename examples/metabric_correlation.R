# examples/metabric_correlation.R

run_example_metabric <- function(input_csv) {
  df <- read.csv(input_csv)
  p <- scatter_plot_with_fit(df, "GeneX", "GeneY")
  ggplot2::ggsave("examples_metabric_corr.png", p)
}
