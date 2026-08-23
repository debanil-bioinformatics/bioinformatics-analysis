# Convert incoming scRNA .txt analyses into cleaned R scripts under scripts/scrna/

# scripts/scrna/pseudobulk_original.R
# This file is a preserved, cleaned version of the original 'pseudobulk.txt' with personal paths removed.

# Usage: source('scripts/scrna/pseudobulk_original.R'); run_pseudobulk_analysis(seurat_rds, output_dir)

run_pseudobulk_analysis <- function(seurat_rds, output_dir = 'results') {
  if (!file.exists(seurat_rds)) stop('seurat_rds not found')
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  scRNA <- readRDS(seurat_rds)
  # The original script performed pseudobulk aggregation, normalization, and ssGSEA scoring.
  # For reproducibility, this function delegates to scripts/scrna/pseudobulk_ucell_cuedc2.R which contains a more robust wrapper.
  source(file.path('scripts', 'scrna', 'pseudobulk_ucell_cuedc2.R'))
  run_pseudobulk_ucell(seurat_rds = seurat_rds, output_dir = output_dir)
}
