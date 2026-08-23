# examples/scrna_signature_scoring.R

# Minimal example showing single-cell signature scoring workflow
run_example_scrna <- function(seurat_rds) {
  seu <- readRDS(seurat_rds)
  seu <- basic_seurat_qc(seu)
  pb <- pseudobulk_aggregate(seu)
  message('Pseudobulk matrix dimensions: ', paste(dim(pb), collapse = ' x '))
}
