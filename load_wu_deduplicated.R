# ==============================================================================
# LOAD WU ET AL. scRNA-seq - FIXED (Remove Duplicates)
# ==============================================================================
# Removes duplicate gene entries before loading
# Memory-optimized batch processing
# ==============================================================================

library(tidyverse)
library(Matrix)
library(Seurat)

cat("\n", strrep("=", 100), "\n")
cat("LOADING WU ET AL. (2021) - DUPLICATE REMOVAL VERSION\n")
cat(strrep("=", 100), "\n\n")

# ==============================================================================
# SETUP
# ==============================================================================

working_dir <- path.expand("~/Documents/Cancers/metadata")
wu_dir <- file.path(working_dir, "GSE161529_Wu_et_al_TNBC")

setwd(wu_dir)

cat("Working directory:", getwd(), "\n\n")

# ==============================================================================
# STEP 1: LOAD FEATURES AND REMOVE DUPLICATES
# ==============================================================================

cat(strrep("-", 100), "\n")
cat("STEP 1: LOADING FEATURES\n")
cat(strrep("-", 100), "\n\n")

features_file <- list.files(pattern = "features.*\\.tsv$", full.names = TRUE, ignore.case = TRUE)[1]

cat("Features file:", features_file, "\n")
features_df <- read.delim(features_file, header = FALSE, stringsAsFactors = FALSE)
colnames(features_df) <- c("gene_id", "gene_name", "feature_type")

cat("Total genes in file:", nrow(features_df), "\n\n")

# Remove duplicates - keep FIRST occurrence of each gene_id
features_df_unique <- features_df[!duplicated(features_df$gene_id), ]

cat("After removing duplicates by gene_id:", nrow(features_df_unique), "\n")
cat("Genes removed:", nrow(features_df) - nrow(features_df_unique), "\n\n")

# Also handle gene_name duplicates by keeping first occurrence
gene_names_orig <- features_df$gene_name
n_dup_names <- sum(duplicated(gene_names_orig))

cat("Duplicate gene NAMES:", n_dup_names, "\n\n")

# Use ENSG IDs as primary identifier (always unique after de-duplication)
gene_ids <- features_df_unique$gene_id

cat("Using", length(gene_ids), "unique ENSG IDs for analysis\n\n")

# ==============================================================================
# STEP 2: FIND FILES
# ==============================================================================

cat(strrep("-", 100), "\n")
cat("STEP 2: FINDING DATA FILES\n")
cat(strrep("-", 100), "\n\n")

matrix_files <- sort(list.files(pattern = "matrix\\.mtx$", full.names = TRUE))
barcode_files <- sort(list.files(pattern = "barcodes\\.tsv$", full.names = TRUE))

cat("Found", length(matrix_files), "samples\n")
sample_names <- sub("_matrix\\.mtx$", "", basename(matrix_files))

cat("Sample names (first 5):\n")
for (i in 1:min(5, length(sample_names))) {
  cat("  ", i, ". ", sample_names[i], "\n", sep="")
}
if (length(sample_names) > 5) {
  cat("  ... and", length(sample_names) - 5, "more\n")
}

cat("\n")

# ==============================================================================
# STEP 3: LOAD SAMPLES IN BATCHES
# ==============================================================================

cat(strrep("-", 100), "\n")
cat("STEP 3: LOADING SAMPLES (BATCH PROCESSING)\n")
cat(strrep("-", 100), "\n\n")

cat("Loading sample 1 as base object...\n")

# Load first sample
mat1 <- Matrix::readMM(matrix_files[1])
bc1 <- read.delim(barcode_files[1], header = FALSE, stringsAsFactors = FALSE)[[1]]

cat("  Matrix dimensions: ", nrow(mat1), " genes x ", ncol(mat1), " cells\n", sep="")

# Subset to only non-duplicate genes
# Keep genes that are in our cleaned gene list
genes_to_keep <- which(1:nrow(mat1) %in% which(!duplicated(features_df$gene_id)))
mat1 <- mat1[genes_to_keep, ]

cat("  After removing duplicates: ", nrow(mat1), " genes x ", ncol(mat1), " cells\n", sep="")

# Set row and column names
rownames(mat1) <- gene_ids
colnames(mat1) <- paste0(sample_names[1], "_", bc1)

# Create initial Seurat object
scRNA <- CreateSeuratObject(
  counts = mat1,
  project = "Wu_et_al_breast",
  min.cells = 3,
  min.features = 200
)

scRNA$sample <- sample_names[1]
scRNA$subtype <- sub("-.*", "", sample_names[1])

cat("  Seurat object created with", ncol(scRNA), "cells\n\n")

# Clean up
rm(mat1, bc1)
gc()

# ==============================================================================
# Load remaining samples and merge incrementally
# ==============================================================================

cat("Adding remaining samples...\n\n")

for (i in 2:length(matrix_files)) {
  
  if (i %% 10 == 0 || i == 2) {
    cat("[", i, "/", length(matrix_files), "] Loading and merging...\n", sep="")
  }
  
  tryCatch({
    
    # Load sample
    mat_i <- Matrix::readMM(matrix_files[i])
    bc_i <- read.delim(barcode_files[i], header = FALSE, stringsAsFactors = FALSE)[[1]]
    
    # Subset to non-duplicate genes
    genes_to_keep <- which(1:nrow(mat_i) %in% which(!duplicated(features_df$gene_id)))
    mat_i <- mat_i[genes_to_keep, ]
    
    # Set names
    rownames(mat_i) <- gene_ids
    colnames(mat_i) <- paste0(sample_names[i], "_", bc_i)
    
    # Create temporary Seurat object
    seurat_i <- CreateSeuratObject(
      counts = mat_i,
      project = "Wu_et_al_breast",
      min.cells = 3,
      min.features = 200
    )
    
    seurat_i$sample <- sample_names[i]
    seurat_i$subtype <- sub("-.*", "", sample_names[i])
    
    # Merge with main object
    scRNA <- merge(scRNA, seurat_i)
    
    # Clean up
    rm(mat_i, bc_i, seurat_i)
    
    if (i %% 10 == 0) {
      gc()
    }
    
  }, error = function(e) {
    cat("  WARNING: Error loading sample", i, "-", e$message, "\n")
  })
}

cat("\n")

# ==============================================================================
# STEP 4: ADD GENE METADATA
# ==============================================================================

cat(strrep("-", 100), "\n")
cat("STEP 4: ADDING GENE METADATA\n")
cat(strrep("-", 100), "\n\n")

# Create mapping
gene_mapping <- setNames(features_df_unique$gene_name, features_df_unique$gene_id)

# Add to object
scRNA@misc$gene_mapping <- gene_mapping

cat("Added gene name mapping\n")
cat("Total genes in object:", nrow(scRNA), "\n\n")

# ==============================================================================
# STEP 5: COMPUTE QC METRICS
# ==============================================================================

cat(strrep("-", 100), "\n")
cat("STEP 5: COMPUTING QC METRICS\n")
cat(strrep("-", 100), "\n\n")

scRNA[["percent.mt"]] <- PercentageFeatureSet(scRNA, pattern = "^MT-|^mt-")
scRNA[["nUMI"]] <- colSums(GetAssayData(scRNA, slot = "counts"))
scRNA[["nGenes"]] <- colSums(GetAssayData(scRNA, slot = "counts") > 0)

cat("QC Metrics Summary:\n\n")
cat("nUMI (library size):\n")
print(summary(scRNA$nUMI))

cat("\nnGenes (genes per cell):\n")
print(summary(scRNA$nGenes))

cat("\nPercent.mt (mitochondrial %):\n")
print(summary(scRNA$percent.mt))

cat("\n")

# ==============================================================================
# STEP 6: SAVE OUTPUT
# ==============================================================================

cat(strrep("-", 100), "\n")
cat("STEP 6: SAVING DATA\n")
cat(strrep("-", 100), "\n\n")

output_dir <- file.path(wu_dir, "processed")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

output_rds <- file.path(output_dir, "Wu_et_al_scRNA_merged.rds")
cat("Saving Seurat object to:", output_rds, "\n")
saveRDS(scRNA, output_rds)

cat("Saving cell metadata...\n")
write.csv(scRNA@meta.data, 
          file.path(output_dir, "cell_metadata.csv"),
          row.names = TRUE)

cat("Saving gene mapping...\n")
write.csv(data.frame(gene_id = names(gene_mapping), 
                     gene_name = unname(gene_mapping)),
          file.path(output_dir, "gene_mapping.csv"),
          row.names = FALSE)

cat("\n")

# ==============================================================================
# SUMMARY
# ==============================================================================

cat(strrep("=", 100), "\n")
cat("✅ LOADING COMPLETE!\n")
cat(strrep("=", 100), "\n\n")

cat("Final Seurat object:\n")
print(scRNA)

cat("\nOutput saved to:", output_dir, "\n\n")

cat("Sample distribution:\n")
print(table(scRNA$sample))

cat("\nSubtype distribution (N, B1, TN, HER2, ER, mER):\n")
print(table(scRNA$subtype))

cat("\nGene information:\n")
cat("  Total genes:", nrow(scRNA), "\n")
cat("  Total cells:", ncol(scRNA), "\n")
cat("  Rownames: ENSG IDs\n")
cat("  Gene symbols available in: scRNA@misc$gene_mapping\n\n")

cat("To access gene names:\n")
cat("  symbol <- scRNA@misc$gene_mapping['ENSG00000163735']  # CUEDC2\n")
cat("  symbol <- scRNA@misc$gene_mapping['ENSG00000091831']  # ESR1\n\n")

cat("To reload in future:\n")
cat("  scRNA <- readRDS('", output_rds, "')\n\n", sep="")
