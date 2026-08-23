# visualization.R - plotting helpers

# simple publication theme
theme_publication <- function(base_size = 12) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
}

# safe labeller for plots
safe_label <- function(text) {
  if (is.null(text)) "" else as.character(text)
}
