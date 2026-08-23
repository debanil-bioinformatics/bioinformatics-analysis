# association.R - clinical association helpers

# compute_contingency_association: chi-square test with simulated p if needed
compute_contingency_association <- function(data, group_col, var_col, group_levels = NULL, simulate.p.value = TRUE) {
  stopifnot(all(c(group_col, var_col) %in% colnames(data)))
  sub <- data[stats::complete.cases(data[[group_col]], data[[var_col]]), , drop = FALSE]
  if (!is.null(group_levels)) sub[[group_col]] <- factor(sub[[group_col]], levels = group_levels)
  tbl <- table(sub[[var_col]], sub[[group_col]])
  if (ncol(tbl) < 2 || nrow(tbl) < 1) return(list(table = tbl, p.value = NA))
  test <- suppressWarnings(chisq.test(tbl, simulate.p.value = simulate.p.value))
  pval <- test$p.value
  perc <- sweep(tbl, 1, rowSums(tbl), "/") * 100
  list(table = tbl, percentages = perc, p.value = pval)
}
