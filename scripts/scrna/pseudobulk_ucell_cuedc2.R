# scripts/scrna/pseudobulk_ucell_cuedc2.R

# Converted from incoming "Pseudobulk and UCell SCoring.txt"; reorganized as a dataset-specific workflow
suppressPackageStartupMessages({
  library(Seurat)
  library(UCell)
  library(dplyr)
})

run_pseudobulk_ucell <- function(seurat_rds, output_dir = "results", gene_sets = NULL) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  seu <- readRDS(seurat_rds)
  gene_map <- seu@misc$gene_mapping
  symbol_to_id <- setNames(names(gene_map), gene_map)
  get_id <- function(symbol) {
    id <- symbol_to_id[symbol]
    id <- id[!is.na(id) & id %in% rownames(seu)]
    if (length(id) == 0) return(NA_character_)
    id[1]
  }
  # prepare gene sets if not provided
  if (is.null(gene_sets)) {
    if (requireNamespace("msigdbr", quietly = TRUE)) {
      hallmark <- msigdbr::msigdbr(species = "Homo sapiens", category = "H")
      early_symbols <- hallmark %>% filter(gs_name == 'HALLMARK_ESTROGEN_RESPONSE_EARLY') %>% pull(gene_symbol) %>% unique()
      late_symbols <- hallmark %>% filter(gs_name == 'HALLMARK_ESTROGEN_RESPONSE_LATE') %>% pull(gene_symbol) %>% unique()
    } else {
      early_symbols <- c('GREB1','PGR','TFF1','ESR1')
      late_symbols <- c('GREB1','ESR1','PGR')
    }
    gene_sets <- list(Estrogen_Response_Early = early_symbols, Estrogen_Response_Late = late_symbols)
  }
  # aggregate pseudobulk
  pseudobulk_counts <- Seurat::AggregateExpression(seu, group.by = 'sample', assays = 'RNA', return.seurat = FALSE)$RNA
  pseudobulk_log <- pseudobulk_counts
  if (requireNamespace("GSVA", quietly = TRUE)) {
    ssgsea <- GSVA::gsva(as.matrix(pseudobulk_log), gene_sets, method = "ssgsea", verbose = FALSE)
    ssgsea <- t(ssgsea)
    write.csv(ssgsea, file.path(output_dir, 'ssgsea_pseudobulk.csv'))
  }
  saveRDS(seu, file.path(output_dir, 'scored_seurat.rds'))
}
