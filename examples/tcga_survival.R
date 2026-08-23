# examples/tcga_survival.R

# Example demonstrating km_fit_and_plot usage
suppressPackageStartupMessages({
  library(readr)
})

example_tcga_survival <- function(example_csv) {
  df <- readr::read_csv(example_csv, show_col_types = FALSE)
  res <- km_fit_and_plot(df, time = "OS_time", event = "OS_event", group = "Biomarker_Status", group_levels = c("Low","High"), max_followup = 200, title = "Example KM")
  ggplot2::ggsave("examples_tcga_km.png", res$plot, width = 6, height = 5)
}
