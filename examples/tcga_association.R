# examples/tcga_association.R

# Example demonstrating compute_contingency_association
suppressPackageStartupMessages({
  library(readr)
})

example_tcga_association <- function(example_csv) {
  df <- readr::read_csv(example_csv, show_col_types = FALSE)
  out <- compute_contingency_association(df, group_col = "Status", var_col = "Age_Group", group_levels = c("Both.High","Others"))
  print(out$p.value)
}
