# Bioinformatics Analysis Workflows

This repository contains bioinformatics analysis scripts developed for the analysis of
large-scale transcriptomic datasets, ChIP-seq experiments, and drug sensitivity assays,
with a focus on cancer biology and hormone receptor signaling.

The code is intended to demonstrate analytical strategies, data handling approaches,
and visualization methods commonly used in computational cancer research.

## Analysis domains
- Bulk RNA-seq analysis and differential gene expression
- Public cohort analysis (TCGA, GEO)
- Survival analysis and biomarker evaluation
- ChIP-seq peak annotation and transcription start site (TSS) proximity analysis
- Large-scale drug sensitivity and IC₅₀ analysis (e.g., GDSC, DepMap)

## Tools and languages
- R (DESeq2, edgeR, limma, survival, GenomicRanges, ggplot2, dplyr)
- Public datasets: TCGA, GEO, GDSC, DepMap

## Repository structure
Each analysis domain is organized into dedicated subdirectories with documented scripts
and minimal assumptions about data structure to support reuse and reproducibility.

## Notes
- Scripts are modular and commented for clarity
- No raw or patient-identifiable data are included
- File paths are user-defined

This repository serves as a living collection of reproducible bioinformatics workflows
rather than a single end-to-end pipeline.
