# Bioinformatics Analysis Toolkit

A lightweight collection of reusable R functions and dataset workflows for transcriptomics, single-cell, and drug-sensitivity analyses. This repository is maintained as a research code portfolio by a computational biologist and contains code used in published and in-progress analyses.

What this repository contains

- R/: reusable analysis functions (survival, correlation, association, expression, signatures, single-cell helpers, visualization)
- scripts/: dataset-specific workflows and wrappers (tcga, metabric, ccle, scrna)
- examples/: short demonstrations showing how to apply the reusable functions to dataset inputs
- tests/: lightweight testthat unit tests for core reusable functions
- docs/: inventory and notes

Design and intent

- R/ holds small, well-documented functions that perform common analytical tasks (e.g., Kaplan–Meier fitting, Cox regression, correlation, signature scoring). Functions accept data frames, matrices or file paths described in the examples; they do not embed dataset-specific assumptions.
- scripts/ holds dataset workflows (TCGA, METABRIC, CCLE/DepMap, Wu et al. scRNA) that call R/ functions and perform the dataset-specific preprocessing needed for those cohorts.
- examples/ contains concise, runnable examples you can adapt to your own data. Examples intentionally use relative paths and accept file paths as arguments.

Quick start

1. Clone the repository and switch to the toolkit/refactor branch:

   git clone https://github.com/debanil-bioinformatics/bioinformatics-analysis.git
   cd bioinformatics-analysis
   git checkout toolkit/refactor

2. Run lightweight tests (recommended):

   # from shell, after installing minimal R dependencies
   Rscript -e "install.packages(c('testthat','survival','survminer','ggplot2'), repos='https://cloud.r-project.org')"
   Rscript -e "testthat::test_dir('tests/testthat')"

3. Try an example (adapt the input file paths):

   Rscript -e "source('examples/tcga_survival.R'); example_tcga_survival('path/to/survival_data.csv')"

Notes

- This repository is research code. It is intentionally lightweight and not a CRAN package. DESCRIPTION is present for convenient metadata only.
- Heavy dependencies used by some scripts (Seurat, GSVA, UCell, DESeq2, biomaRt) are optional — install them only if you plan to run the corresponding scripts in scripts/scrna/ or scripts/tcga/.
- Original raw analysis scripts and prior versions are available in the repository history. Current, canonical workflow scripts are in scripts/ and reusable functions are in R/.

If you find any issues running the examples or tests, please report them so they can be addressed; the repository is intended to be useful and reproducible for other computational biology researchers.
