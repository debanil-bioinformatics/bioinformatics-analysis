############################################################
# Author: Debanil Dhar
# Year: 2025
#
# Description:
# This script downloads raw gene-level RNA-seq counts
# (STAR - Counts workflow) from multiple TCGA projects
# using the TCGAbiolinks package and saves each project
# as an RDS file for downstream transcriptomic analysis.
#
# Input:
# - TCGA project identifiers
#
# Output:
# - One RDS file per TCGA project containing raw count data
#
# Notes:
# - No patient-identifiable data are included
# - Designed for batch download and reproducible reuse
############################################################

## Load required libraries
suppressPackageStartupMessages({
  library(TCGAbiolinks)
})

## Define output directory
output_dir <- "tcga_star_counts"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

## List of TCGA projects (excluding TCGA-UCEC)
cancer_projects <- c(
  "TCGA-BLCA", "TCGA-BRCA", "TCGA-CESC", "TCGA-CHOL", "TCGA-COAD",
  "TCGA-ESCA", "TCGA-GBM", "TCGA-HNSC", "TCGA-KICH", "TCGA-KIRC",
  "TCGA-KIRP", "TCGA-LIHC", "TCGA-LUAD", "TCGA-LUSC", "TCGA-PAAD",
  "TCGA-PRAD", "TCGA-PCPG", "TCGA-READ", "TCGA-SARC", "TCGA-SKCM",
  "TCGA-THCA", "TCGA-THYM", "TCGA-STAD"
)

## Loop through each project and download data
for (project in cancer_projects) {
  
  message("Processing project: ", project)
  
  query <- GDCquery(
    project = project,
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts"
  )
  
  GDCdownload(query)
  data <- GDCprepare(query)
  
  output_file <- file.path(
    output_dir,
    paste0(project, "_STAR_counts.rds")
  )
  
  saveRDS(data, file = output_file)
  message("Saved: ", output_file)
}

message("All TCGA projects downloaded successfully.")
