# Fix examples/ccle_correlation.R to be a concise functional example

# Example: run_ccle_bub1b_mad2l1
# Usage: source('examples/ccle_correlation.R'); run_example_ccle(expr_file, meta_file, output_dir)

run_example_ccle <- function(expr_file, meta_file, output_dir = 'results') {
  if (!file.exists(expr_file) || !file.exists(meta_file)) stop('Input files not found. Provide expression CSV and metadata CSV paths.')
  scripts::ccle <- NULL
  # call the wrapper implemented in scripts/ccle/correlation_ccle.R
  source(file.path('scripts', 'ccle', 'correlation_ccle.R'))
  run_ccle_bub1b_mad2l1(expr_file = expr_file, meta_file = meta_file, output_dir = output_dir)
}
