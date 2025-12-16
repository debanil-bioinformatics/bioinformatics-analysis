# Transcriptomics Analysis

This directory contains reproducible bioinformatics workflows for bulk RNA-seq analysis
using large public cancer cohorts, with an emphasis on differential expression,
clinical integration, and survival analysis.

## Analysis scope
- **TCGA data acquisition and preprocessing**  
  Batch download and handling of RNA-seq STAR-count data using TCGAbiolinks.

- **Differential gene expression analysis**  
  Gene-level differential expression using DESeq2, including gene subset–based analyses.

- **Cross-cohort gene ranking**  
  Integration and rank aggregation of differential expression results across multiple
  TCGA cancer types.

- **Clinicopathological association analysis**  
  Association testing between molecular features and clinical variables using
  contingency analysis and publication-ready tabulation.

- **Multivariate Cox proportional hazards modeling**  
  Adjusted survival analysis incorporating clinical covariates across TCGA and
  METABRIC cohorts.

- **Kaplan–Meier survival analysis**  
  Generalized survival curve generation with hazard ratio estimation and
  cohort-agnostic visualization.

## Notes
- Scripts are designed to be **dataset-agnostic and reusable**
- No raw or patient-identifiable data are included
- All analyses are implemented in R using established Bioconductor and CRAN packages
