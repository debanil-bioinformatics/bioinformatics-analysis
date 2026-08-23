# Final cleanup commit: delete obsolete pointer files and tidy structure

This commit removes legacy pointer files that were replaced by canonical scripts under scripts/ and reusable functions in R/. It also deletes obsolete .txt copies of R scripts whose content was migrated.

Files removed:
- TCGA Correlation.R, CCLE Correlation.R, METABRIC Correlation.R, TCGA Regression.R,
  48 months KM.R, HR- TCGA Association.R, Tumor vs Normal.R, UV MV.R, tcga_association_1.R,
  load_wu_deduplicated.R, metabric_association.R, metabric_cox_complete_3.R,
  Pseudobulk and UCell SCoring.txt, pseudobulk.txt, scRNA final analysis.txt,
  scRNA seq clustering Wu Et Al.txt, scRNA step 6 QC.txt

These deletions are safe because the canonical code is present under scripts/ and R/.
