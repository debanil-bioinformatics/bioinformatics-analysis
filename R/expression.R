# expression.R - expression matrix helpers

# remove_ensembl_version: drop version suffix from Ensembl IDs
remove_ensembl_version <- function(ids) {
  sub("\\..*$", "", ids)
}

# subset_genes: subset matrix or dataframe by gene ids or symbols
subset_genes <- function(expr, genes) {
  stopifnot(is.matrix(expr) || is.data.frame(expr))
  genes_present <- intersect(rownames(expr), genes)
  expr[genes_present, , drop = FALSE]
}

# map_ensembl_to_symbol: wrapper around biomaRt for mapping; returns named vector
map_ensembl_to_symbol <- function(ensembl_ids, host = "https://www.ensembl.org") {
  if (!requireNamespace("biomaRt", quietly = TRUE)) stop("biomaRt required for mapping")
  mart <- biomaRt::useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", host = host)
  res <- biomaRt::getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                        filters = "ensembl_gene_id", values = ensembl_ids, mart = mart)
  sym <- res$hgnc_symbol
  names(sym) <- res$ensembl_gene_id
  sym
}
