# single_cell.R - Seurat/SCRNA helpers

# basic_seurat_qc: perform common QC filters and normalization
basic_seurat_qc <- function(seu, min.features = 200, max.features = NULL, min.cells = 3, normalize.method = "LogNormalize") {
  stopifnot(requireNamespace("Seurat", quietly = TRUE))
  # filter cells
  if (!is.null(min.features)) seu <- Seurat::subset(seu, subset = nFeature_RNA >= min.features)
  if (!is.null(max.features)) seu <- Seurat::subset(seu, subset = nFeature_RNA <= max.features)
  # filter genes
  seu <- Seurat::subset(seu, features = rownames(seu)[Matrix::rowSums(Seurat::GetAssayData(seu, assay = "RNA", slot = "counts") > 0) >= min.cells])
  seu <- Seurat::NormalizeData(seu, normalization.method = normalize.method)
  seu <- Seurat::FindVariableFeatures(seu)
  seu
}

# pseudobulk_aggregate: aggregate counts by sample metadata column
pseudobulk_aggregate <- function(seu, sample_col = "sample", assay = "RNA", slot = "counts") {
  stopifnot(requireNamespace("Seurat", quietly = TRUE))
  mat <- Seurat::AggregateExpression(seu, group.by = sample_col, assays = assay, slot = slot, return.seurat = FALSE)
  # returns list with assay name
  mat[[assay]]
}
