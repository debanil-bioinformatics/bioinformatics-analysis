# correlation.R - reusable correlation helpers

# gene_correlation: compute correlation between two numeric vectors (by name)
gene_correlation <- function(x, y, method = "spearman", use = "complete.obs") {
  stopifnot(is.numeric(x), is.numeric(y))
  res <- stats::cor.test(x, y, method = method, exact = FALSE)
  list(r = unname(res$estimate), p.value = res$p.value, method = method)
}

# correlation_table: apply correlation across pairs and return tibble
correlation_table <- function(df, x, y, method = "spearman") {
  stopifnot(all(c(x, y) %in% colnames(df)))
  xvec <- df[[x]]
  yvec <- df[[y]]
  r <- gene_correlation(xvec, yvec, method = method)
  data.frame(x = x, y = y, r = r$r, p.value = r$p.value, n = sum(stats::complete.cases(xvec, yvec)))
}

# scatter_plot_with_fit: simple scatter with linear fit and annotation
scatter_plot_with_fit <- function(df, x, y, xlab = NULL, ylab = NULL, title = NULL, method = "pearson") {
  stopifnot(all(c(x, y) %in% colnames(df)))
  res <- gene_correlation(df[[x]], df[[y]], method = method)
  label <- sprintf("r = %.3f\np = %.2e\nn = %d", res$r, res$p.value, sum(stats::complete.cases(df[[x]], df[[y]])))
  p <- ggplot2::ggplot(df, ggplot2::aes_string(x = x, y = y)) +
    ggplot2::geom_point(alpha = 0.6) +
    ggplot2::geom_smooth(method = "lm", se = TRUE) +
    ggplot2::annotate("text", x = Inf, y = -Inf, label = label, hjust = 1.1, vjust = -0.5, size = 3, family = "mono") +
    ggplot2::labs(title = title, x = xlab %||% x, y = ylab %||% y) +
    ggplot2::theme_classic()
  p
}
