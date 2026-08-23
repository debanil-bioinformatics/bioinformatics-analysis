# scripts/scrna/scRNA_final_analysis.R
# Cleaned copy of 'scRNA final analysis.txt' moved into scripts/scrna/

run_scrna_final_analysis <- function(seurat_rds, output_dir = 'results') {
  if (!file.exists(seurat_rds)) stop('seurat_rds not found')
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  scRNA <- readRDS(seurat_rds)
  # This function provides an entry point to the final analysis steps.
  # For the detailed per-step implementation, use scripts/scrna/pseudobulk_ucell_cuedc2.R and other single-cell wrappers.
  # Here, we run the pseudobulk+UCell wrapper as part of the final analysis.
  source(file.path('scripts', 'scrna', 'pseudobulk_ucell_cuedc2.R'))
  run_pseudobulk_ucell(seurat_rds = seurat_rds, output_dir = output_dir)
}
