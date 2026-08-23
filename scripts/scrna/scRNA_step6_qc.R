# scripts/scrna/scRNA_step6_qc.R
# Cleaned copy of 'scRNA step 6 QC.txt' moved into scripts/scrna/

run_scrna_qc <- function(seurat_rds, output_dir = 'results', min.features = 200, max.features = NULL) {
  if (!file.exists(seurat_rds)) stop('seurat_rds not found')
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  seu <- readRDS(seurat_rds)
  if (!requireNamespace('Seurat', quietly = TRUE)) stop('Seurat required for QC')
  seu <- basic_seurat_qc(seu, min.features = min.features, max.features = max.features)
  saveRDS(seu, file.path(output_dir, 'seu_qc.rds'))
}
