# Repository inventory for refactor

This file is an inventory created to support the toolkit/refactor work. I inspected the top-level tree and a selection of relevant source files (R / .txt) to build this summary. No files were modified.

## Top-level structure (current)

- .gitignore
- LICENSE
- README.md
- Several root-level R scripts and R-like text files (see list below)
- chipseq/ (README.md)
- drug_sensitivity/ (README.md)
- transcriptomics/
  - README.md
  - generalized_kaplan_meier_and_cox.R
  - generalized_multivariate_cox_regression.R
  - tcga_clinicopathological_association.R
  - tcga_gene_subset_deseq2.R
  - tcga_rank_aggregation_across_projects.R
  - tcga_star_counts_download.R

### Newly uploaded / incoming root files (R and .txt)
- 48 months KM.R
- CCLE Correlation.R
- `HR- TCGA Association.R`
- METABRIC Correlation.R
- `Pseudobulk and UCell SCoring.txt`
- TCGA Correlation.R
- TCGA Regression.R
- Tumor vs Normal.R
- UV MV.R
- load_wu_deduplicated.R
- metabric_association.R
- metabric_cox_complete_3.R
- pseudobulk.txt
- scRNA final analysis.txt
- scRNA seq clustering Wu Et Al.txt
- scRNA step 6 QC.txt
- tcga_association_1.R

(These are all at repository root in the current main branch.)

## Existing analysis areas

- Transcriptomics (existing `transcriptomics/` folder)
  - Survival analysis (KM, Cox) scripts
  - Multivariate Cox regression with gtsummary export
  - Clinicopathological association (chi-square) for TCGA
  - DESeq2 gene-subset analysis for TCGA
  - STAR counts handling / download helper
- ChIP-seq (chipseq/) — only README present
- Drug sensitivity (drug_sensitivity/) — only README present
- Single-cell (incoming files) — multiple Wu et al. scRNA-seq analysis scripts in .txt files
- Cross-dataset correlation/regression workflows (root-level correlation scripts for TCGA / METABRIC / CCLE)

## Reusable functions already present (identified by name and file)

- transcriptomics/tcga_clinicopathological_association.R
  - compute_association(data, var, display_name) — computes contingency table, runs (simulated) chi-square, returns p-value and formatted results.

- `Pseudobulk and UCell SCoring.txt` (incoming) contains several small helper functions used in-scene:
  - get_id(symbol) — map gene symbol to ENSG from a Seurat object gene_map
  - run_ucell_corr(df, group_label) — perform Spearman correlations on UCell scores and return tibble
  - run_ucell_wilcox(df, group_label) — median-split Wilcoxon comparisons and return tibble
  (These are implemented as inline helpers in the analysis script and are prime candidates for extraction.)

(No dedicated R package-style utility files (R/) are present yet.)

## Dataset-specific scripts (candidate workflows to move to `scripts/`)

- transcriptomics/tcga_* (multiple): clearly TCGA dataset workflows
- root-level TCGA Correlation.R, TCGA Regression.R, tcga_association_1.R, HR- TCGA Association.R: TCGA-specific analyses
- CCLE Correlation.R: CCLE/DepMap correlation workflow (cell-line filtering and plotting)
- METABRIC Correlation.R, metabric_association.R, metabric_cox_complete_3.R: METABRIC-specific analyses
- Several scRNA-related .txt files (Wu et al. TNBC workflows): pseudobulk, UCell scoring, QC, clustering — dataset-specific scRNA workflows
- load_wu_deduplicated.R: loader / helper for Wu dataset
- `48 months KM.R`, `Tumor vs Normal.R`, `UV MV.R`, `UV MV.R` (likely univariable/multivariable) — small root-level scripts that perform survival or group-based comparisons (require inspection)

## `.txt` files that contain R code (must be treated as source)
- Pseudobulk and UCell SCoring.txt
- pseudobulk.txt
- scRNA final analysis.txt
- scRNA seq clustering Wu Et Al.txt
- scRNA step 6 QC.txt

These appear to be full R scripts saved with .txt extension and should be inspected and renamed to `.R` if incorporated.

## Obvious duplicates / potentially obsolete files

- Multiple TCGA association/regression scripts at root and in transcriptomics/ (tcga_clinicopathological_association.R, tcga_association_1.R, `HR- TCGA Association.R`, TCGA Regression.R). These likely overlap in functionality and will need consolidation.
- Multiple correlation scripts for different datasets (TCGA Correlation.R, CCLE Correlation.R, METABRIC Correlation.R) implement very similar correlation+plotting logic; candidate for consolidation into a reusable correlation function + small dataset-specific wrapper scripts.
- Transcriptomics contains two generalized survival scripts (`generalized_kaplan_meier_and_cox.R`, `generalized_multivariate_cox_regression.R`) plus several root-level KM/regression scripts (e.g., `48 months KM.R`, `UV MV.R`, `TCGA Regression.R`) — consolidate shared logic.

I did NOT delete or modify any files; this is only an inventory based on file names and targeted content inspection.

## Major repeated analytical operations (good targets for R/ functions)

- Correlation (Pearson/Spearman), including tidy output of r, p, R², n and publication-quality scatter plots.
- Kaplan–Meier fitting and plotting + Cox regression (univariable and multivariable), extracting HR/p, risk tables, and publication-ready plotting.
- Clinical association (contingency tables / chi-square with simulated p-values and formatted output tables).
- Differential expression (DESeq2) for gene subsets and general expression matrix preprocessing (ENSEMBL mapping, version-stripping).
- Ensembl↔HGNC mapping (biomaRt usage) used in DESeq2 script.
- Single-cell reusable steps: loading Seurat objects, QC filtering, pseudobulk aggregation, normalization, ssGSEA/GSVA and UCell scoring, chunked UCell processing, pseudobulk scoring, per-sample aggregation and tests (Spearman, Wilcoxon), nUMI/nGene confound checks, plotting with label repulsion.
- Signature scoring (ssGSEA via GSVA; UCell; simple mean/median scoring) across bulk and single-cell data.
- TCGA-specific barcode parsing (sample type inference from barcode parts) — dataset-specific helper logic.

## Files that should probably be preserved unchanged (initially)

- LICENSE (MIT)
- README.md (will be rewritten, but preserve content/history)
- transcriptomics/*.R core scripts until their logic is extracted into R/ (don't delete until new functions exist and examples reproduce the behavior)
- scRNA incoming analysis .txt files — preserve content; they contain useful, tested analysis flows (Wu et al.) and should be migrated into `scripts/scrna/` and examples after extraction of reusable functions.
- load_wu_deduplicated.R — likely dataset loader; preserve as a script or move to scripts/scrna/

## Immediate next steps (not executed yet)

- Create R/ with the targeted reusable files and extract functions from the scripts listed above.
- Move dataset-specific scripts into `scripts/` (tcga/, metabric/, ccle/, scrna/) as lightweight wrappers that call R/ functions.
- Convert `.txt` files containing R code to `.R` and relocate into `scripts/scrna/` or examples/ as appropriate.
- Add `examples/` demonstrating the common workflows (TCGA survival, association, CCLE correlation, scrna signature scoring).
- Add lightweight tests for survival, correlation and signature-scoring functions.
- Update README, create DESCRIPTION, and create docs/ with this inventory and any short how-to.

---

Inventory created programmatically and saved to `docs/refactor_inventory.md` on branch `toolkit/refactor`.

If you'd like, I can now proceed to implement the next changes (extract R/ functions, move scripts, rename .txt→.R, add examples and tests).