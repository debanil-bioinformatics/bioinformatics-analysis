# Updated correlation helpers with safer labs

# gene_correlation: compute correlation between two numeric vectors (by name)
gene_correlation <- function(x, y, method = "spearman", use = "complete.obs") {
  if (!is.numeric(x) || !is.numeric(y)) stop("x and y must be numeric vectors")
  if (length(x) != length(y)) stop("x and y must have the same length")
  res <- stats::cor.test(x, y, method = method, exact = FALSE)
  list(r = unname(res$estimate), p.value = res$p.value, method = method)
}

# correlation_table: apply correlation across pairs and return data.frame
correlation_table <- function(df, x, y, method = "spearman") {
  stopifnot(all(c(x, y) %in% colnames(df)))
  xvec <- as.numeric(df[[x]])
  yvec <- as.numeric(df[[y]])
  r <- gene_correlation(xvec, yvec, method = method)
  data.frame(x = x, y = y, r = r$r, p.value = r$p.value, n = sum(stats::complete.cases(xvec, yvec)), stringsAsFactors = FALSE)
}

# scatter_plot_with_fit: simple scatter with linear fit and annotation
scatter_plot_with_fit <- function(df, x, y, xlab = NULL, ylab = NULL, title = NULL, method = "pearson") {
  stopifnot(all(c(x, y) %in% colnames(df)))
  xvec <- as.numeric(df[[x]])
  yvec <- as.numeric(df[[y]])
  res <- gene_correlation(xvec, yvec, method = method)
  label <- sprintf("r = %.3f\np = %.2e\nn = %d", res$r, res$p.value, sum(stats::complete.cases(xvec, yvec)))
  xlab_use <- if (!is.null(xlab)) xlab else x
  ylab_use <- if (!is.null(ylab)) ylab else y
  p <- ggplot2::ggplot(df, ggplot2::aes_string(x = x, y = y)) +
    ggplot2::geom_point(alpha = 0.6) +
    ggplot2::geom_smooth(method = "lm", se = TRUE) +
    ggplot2::annotate("text", x = Inf, y = -Inf, label = label, hjust = 1.05, vjust = -0.5, size = 3, family = "mono") +
    ggplot2::labs(title = title, x = xlab_use, y = ylab_use) +
    ggplot2::theme_classic()
  p
}
