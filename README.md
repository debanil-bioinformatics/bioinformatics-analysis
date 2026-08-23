# Bioinformatics Analysis Toolkit

A lightweight collection of reusable R functions and dataset workflows for transcriptomics, single-cell, and drug-sensitivity analyses.

Structure

- R/: reusable functions (survival, correlation, association, expression, signatures, single_cell, visualization)
- scripts/: dataset-specific workflows (tcga, metabric, ccle, scrna)
- examples/: concise example scripts
- tests/: lightweight testthat tests

See docs/refactor_inventory.md for a detailed inventory of repository contents collected before refactor.

Dependencies

Key R packages used by functions: survival, survminer, ggplot2, dplyr, DESeq2 (for DE workflows), Seurat (for single-cell), GSVA/UCell (for signature scoring).

Usage

Source the function files in R/ or use devtools::load_all('.') in an interactive session. Run dataset workflows under scripts/ with appropriate input files.
