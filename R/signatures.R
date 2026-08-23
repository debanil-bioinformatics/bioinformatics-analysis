# signatures.R - signature scoring utilities

# score_by_mean: simple signature score = mean expression across genes
score_by_mean <- function(expr_matrix, gene_set) {
  common <- intersect(rownames(expr_matrix), gene_set)
  if (length(common) == 0) return(rep(NA_real_, ncol(expr_matrix)))
  colMeans(expr_matrix[common, , drop = FALSE], na.rm = TRUE)
}

# score_ssgsea: wrapper for GSVA::gsva with method ssgsea
score_ssgsea <- function(expr_matrix, gene_sets) {
  if (!requireNamespace("GSVA", quietly = TRUE)) stop("GSVA required for ssgsea scoring")
  res <- GSVA::gsva(as.matrix(expr_matrix), gene_sets, method = "ssgsea", verbose = FALSE)
  t(res)
}

# score_ucell: wrapper placeholder for UCell scoring (per-cell)
score_ucell <- function(expr_matrix, gene_sets, ...) {
  if (!requireNamespace("UCell", quietly = TRUE)) stop("UCell required for UCell scoring")
  UCell::ScoreSignatures_UCell(expr_matrix, features = gene_sets, ...)
}
