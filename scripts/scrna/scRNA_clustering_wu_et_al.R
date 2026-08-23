# scripts/scrna/scRNA_clustering_wu_et_al.R
# Cleaned copy of 'scRNA seq clustering Wu Et Al.txt' moved into scripts/scrna/

run_scrna_clustering <- function(seurat_rds, output_dir = 'results') {
  if (!file.exists(seurat_rds)) stop('seurat_rds not found')
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  seu <- readRDS(seurat_rds)
  # perform standard Seurat workflow steps using reusable functions where appropriate
  if (!requireNamespace('Seurat', quietly = TRUE)) stop('Seurat required for clustering')
  seu <- Seurat::NormalizeData(seu)
  seu <- Seurat::FindVariableFeatures(seu)
  seu <- Seurat::ScaleData(seu)
  seu <- Seurat::RunPCA(seu)
  seu <- Seurat::FindNeighbors(seu, dims = 1:20)
  seu <- Seurat::FindClusters(seu, resolution = 0.5)
  seu <- Seurat::RunUMAP(seu, dims = 1:20)
  saveRDS(seu, file.path(output_dir, 'seu_clustered.rds'))
}
