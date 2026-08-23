library(survival)
library(broom)
library(dplyr)
library(tidyr)
library(flextable)
library(officer)
library(png)
library(grid)
library(readr)

setwd("~/Documents/Cancers")

cat("Starting METABRIC Cox Regression Analysis...\n\n")

# ── 1. LOAD DATA ──────────────────────────────────────────────────────────────

cat("Loading data...\n")

# Expression data
expr <- read_tsv("brca_metabric/data_mrna_illumina_microarray.txt",
                 show_col_types = FALSE)

# Clinical data
clinical_patient <- read_tsv("brca_metabric/data_clinical_patient.txt",
                             show_col_types = FALSE)
clinical_sample <- read_tsv("brca_metabric/data_clinical_sample.txt",
                            show_col_types = FALSE)

# ── 2. EXTRACT BUB1B AND MAD2L1 EXPRESSION ────────────────────────────────────

bub1b_expr <- expr %>%
  filter(Hugo_Symbol == "BUB1B") %>%
  select(-Hugo_Symbol, -Entrez_Gene_Id) %>%
  pivot_longer(everything(), names_to = "SampleID", values_to = "BUB1B_expr")

mad2l1_expr <- expr %>%
  filter(Hugo_Symbol == "MAD2L1") %>%
  select(-Hugo_Symbol, -Entrez_Gene_Id) %>%
  pivot_longer(everything(), names_to = "SampleID", values_to = "MAD2L1_expr")

# ── 3. PREPARE CLINICAL DATA ──────────────────────────────────────────────────

colnames(clinical_patient) <- gsub(" ", "_", colnames(clinical_patient))
colnames(clinical_sample) <- gsub(" ", "_", colnames(clinical_sample))

clinical <- clinical_sample %>%
  select(`#Patient_Identifier`, 
         Sample_Identifier,
         Tumor_Stage,
         Neoplasm_Histologic_Grade,
         Tumor_Size) %>%
  rename(PatientID = `#Patient_Identifier`,
         SampleID = Sample_Identifier,
         Stage = Tumor_Stage,
         Grade = Neoplasm_Histologic_Grade) %>%
  left_join(
    clinical_patient %>%
      select(`#Patient_Identifier`,
             `Age_at_Diagnosis`,
             `Lymph_nodes_examined_positive`,
             `Inferred_Menopausal_State`,
             `Pam50_+_Claudin-low_subtype`,
             `Overall_Survival_(Months)`,
             `Overall_Survival_Status`,
             `Relapse_Free_Status_(Months)`,
             `Relapse_Free_Status`) %>%
      rename(PatientID = `#Patient_Identifier`,
             Age = `Age_at_Diagnosis`,
             Lymph_Nodes = `Lymph_nodes_examined_positive`,
             Menopausal = `Inferred_Menopausal_State`,
             PAM50 = `Pam50_+_Claudin-low_subtype`,
             OS_Months = `Overall_Survival_(Months)`,
             OS_Status = `Overall_Survival_Status`,
             RFS_Months = `Relapse_Free_Status_(Months)`,
             RFS_Status = `Relapse_Free_Status`),
    by = "PatientID"
  )

# ── 4. CATEGORIZE VARIABLES ───────────────────────────────────────────────────

clinical <- clinical %>%
  mutate(
    Age_Group = factor(
      ifelse(Age < 60, "Low", "High"),
      levels = c("Low", "High")
    ),
    
    Nodal_Status = factor(
      ifelse(Lymph_Nodes <= 1, "Low", "High"),
      levels = c("Low", "High")
    ),
    
    Grade_Binary = factor(
      ifelse(Grade %in% c("1", "2"), "Low", "High"),
      levels = c("Low", "High")
    ),
    
    Stage_Binary = factor(
      ifelse(Stage %in% c("1", "2"), "Low", "High"),
      levels = c("Low", "High")
    ),
    
    Menopausal = factor(Menopausal, levels = c("Pre", "Post")),
    
    # PAM50: Luminal (LumA + LumB) vs Non-Luminal (Her2 + Basal + claudin-low + Normal)
    PAM50_Binary = factor(
      case_when(
        PAM50 %in% c("LumA", "LumB") ~ "Luminal",
        PAM50 %in% c("Her2", "Basal", "claudin-low", "Normal") ~ "Non-Luminal",
        TRUE ~ NA_character_
      ),
      levels = c("Luminal", "Non-Luminal")
    )
  ) %>%
  select(-Age, -Lymph_Nodes, -Tumor_Size, -PAM50, -Grade, -Stage)

# ── 5. MERGE EXPRESSION WITH CLINICAL DATA ────────────────────────────────────

data <- bub1b_expr %>%
  left_join(mad2l1_expr, by = "SampleID") %>%
  left_join(clinical, by = "SampleID")

cat("Total samples before filtering:", nrow(data), "\n")

# ── 6. CONVERT SURVIVAL MONTHS TO NUMERIC ─────────────────────────────────────

data <- data %>%
  mutate(
    OS_Months = as.numeric(OS_Months),
    RFS_Months = as.numeric(RFS_Months)
  )

# ── 7. REMOVE NAS AND CREATE FINAL DATASET ────────────────────────────────────

data <- data %>%
  drop_na(BUB1B_expr, MAD2L1_expr, OS_Months, OS_Status, RFS_Months, RFS_Status,
          Age_Group, Stage_Binary, Nodal_Status, Grade_Binary, Menopausal, PAM50_Binary)

cat("Total samples with complete data:", nrow(data), "\n\n")

# ── 8. CONVERT SURVIVAL STATUS TO NUMERIC ─────────────────────────────────────

cat("Survival Status Values:\n")
cat("  OS_Status unique:", paste(unique(data$OS_Status), collapse = ", "), "\n")
cat("  RFS_Status unique:", paste(unique(data$RFS_Status), collapse = ", "), "\n\n")

data <- data %>%
  mutate(
    OS_Status_numeric = as.numeric(grepl("1|DECEASED|Dead", OS_Status, ignore.case = TRUE)),
    RFS_Status_numeric = as.numeric(grepl("1|Recurred|Relapsed", RFS_Status, ignore.case = TRUE))
  )

# ── 9. CATEGORIZE BUB1B AND MAD2L1 BY MEDIAN ──────────────────────────────────

data <- data %>%
  mutate(
    BUB1B_Status = factor(
      ifelse(BUB1B_expr >= median(BUB1B_expr, na.rm = TRUE), "High", "Low"),
      levels = c("Low", "High")
    ),
    MAD2L1_Status = factor(
      ifelse(MAD2L1_expr >= median(MAD2L1_expr, na.rm = TRUE), "High", "Low"),
      levels = c("Low", "High")
    )
  )

# ── 10. FIT MULTIVARIATE COX MODELS ───────────────────────────────────────────

cat("Fitting Cox proportional hazards models...\n\n")

# Overall Survival Model
cox_os <- coxph(
  Surv(OS_Months, OS_Status_numeric) ~ 
    Age_Group + Stage_Binary + Nodal_Status + Grade_Binary + 
    Menopausal + PAM50_Binary + BUB1B_Status + MAD2L1_Status,
  data = data
)

# Recurrence-Free Survival Model
cox_rfs <- coxph(
  Surv(RFS_Months, RFS_Status_numeric) ~ 
    Age_Group + Stage_Binary + Nodal_Status + Grade_Binary + 
    Menopausal + PAM50_Binary + BUB1B_Status + MAD2L1_Status,
  data = data
)

# ── 11. EXTRACT AND TIDY RESULTS ──────────────────────────────────────────────

# Tidy OS results
tidy_os <- broom::tidy(cox_os, exponentiate = TRUE, conf.int = TRUE) %>%
  mutate(
    variable = term,
    HR_CI = sprintf("%.2f (%.2f-%.2f)", estimate, conf.low, conf.high),
    pvalue_fmt = case_when(
      p.value < 0.001 ~ "<0.001",
      TRUE ~ sprintf("%.3f", p.value)
    )
  ) %>%
  select(variable, HR_CI_OS = HR_CI, pvalue_OS = pvalue_fmt)

# Tidy RFS results
tidy_rfs <- broom::tidy(cox_rfs, exponentiate = TRUE, conf.int = TRUE) %>%
  mutate(
    variable = term,
    HR_CI = sprintf("%.2f (%.2f-%.2f)", estimate, conf.low, conf.high),
    pvalue_fmt = case_when(
      p.value < 0.001 ~ "<0.001",
      TRUE ~ sprintf("%.3f", p.value)
    )
  ) %>%
  select(variable, HR_CI_RFS = HR_CI, pvalue_RFS = pvalue_fmt)

# Merge OS and RFS results
merged_tbl <- full_join(tidy_os, tidy_rfs, by = "variable")

# ── 12. BUILD CLEAN FORMATTED RESULTS TABLE ───────────────────────────────────

var_mapping <- tribble(
  ~term, ~group, ~ref, ~comp,
  "Age_GroupHigh", "Age Group", "Low (<60)", "High (≥60)",
  "Stage_BinaryHigh", "Tumor Stage", "Low (I-II)", "High (III-IV)",
  "Nodal_StatusHigh", "Nodal Involvement", "Low (N0-N1)", "High (N2-N3)",
  "Grade_BinaryHigh", "Histologic Grade", "Low (1-2)", "High (3)",
  "MenopausalPost", "Menopausal State", "Pre", "Post",
  "PAM50_BinaryNon-Luminal", "PAM50 Subtype", "Luminal", "Non-Luminal",
  "BUB1B_StatusHigh", "BUB1B Expression", "Low", "High",
  "MAD2L1_StatusHigh", "MAD2L1 Expression", "Low", "High"
)

# Create clean table structure
result_list <- list()

for (i in 1:nrow(var_mapping)) {
  row_info <- var_mapping[i, ]
  term <- row_info$term
  main_row <- merged_tbl %>% filter(variable == term)
  
  if (nrow(main_row) > 0) {
    # Reference category
    ref_list <- list(
      Variables = row_info$group,
      Category = row_info$ref,
      OS_HR = "1.00",
      OS_P = "—",
      RFS_HR = "1.00",
      RFS_P = "—"
    )
    result_list[[paste0(i, "_ref")]] <- as.data.frame(ref_list, stringsAsFactors = FALSE)
    
    # Comparative category
    comp_list <- list(
      Variables = "",
      Category = paste0("  ", row_info$comp),
      OS_HR = main_row$HR_CI_OS[1],
      OS_P = main_row$pvalue_OS[1],
      RFS_HR = main_row$HR_CI_RFS[1],
      RFS_P = main_row$pvalue_RFS[1]
    )
    result_list[[paste0(i, "_comp")]] <- as.data.frame(comp_list, stringsAsFactors = FALSE)
  }
}

merged_tbl_formatted <- do.call(rbind, result_list) %>%
  as_tibble() %>%
  rename(
    `Overall Survival\nHazard Ratio (95% CI)` = OS_HR,
    `Overall Survival\np-Value` = OS_P,
    `RFS\nHazard Ratio (95% CI)` = RFS_HR,
    `RFS\np-Value` = RFS_P
  )

# ── 13. CREATE PUBLICATION-QUALITY FLEXTABLE ──────────────────────────────────

ft <- flextable(merged_tbl_formatted) %>%
  bold(part = "header") %>%
  bold(i = ~Variables != "", j = "Variables") %>%
  align(align = "left", j = c("Variables", "Category"), part = "body") %>%
  align(align = "center", j = c("Overall Survival\nHazard Ratio (95% CI)", 
                                 "Overall Survival\np-Value",
                                 "RFS\nHazard Ratio (95% CI)",
                                 "RFS\np-Value"), part = "body") %>%
  border_remove() %>%
  hline_top(part = "header", border = officer::fp_border(color = "black", width = 2)) %>%
  hline_bottom(part = "body", border = officer::fp_border(color = "black", width = 2)) %>%
  autofit() %>%
  fit_to_width(11)

# ── 14. SAVE AS PDF ───────────────────────────────────────────────────────────

cat("Saving results to PDF...\n\n")

# Save as PNG first
save_as_image(ft, path = "metabric_cox_results.png", zoom = 2)

# Create PDF
pdf("METABRIC_Multivariate_Cox_Regression_BUB1B_MAD2L1.pdf", width = 12, height = 8)
grid.newpage()
img <- readPNG("metabric_cox_results.png")
grid.raster(img)
dev.off()

# ── 15. PRINT SUMMARY OUTPUT ──────────────────────────────────────────────────

cat("════════════════════════════════════════════════════════════════════════════\n")
cat("✅ MULTIVARIATE COX REGRESSION ANALYSIS COMPLETE!\n")
cat("════════════════════════════════════════════════════════════════════════════\n\n")

cat("PDF saved: METABRIC_Multivariate_Cox_Regression_BUB1B_MAD2L1.pdf\n\n")

cat("─── OVERALL SURVIVAL MODEL ───\n")
print(summary(cox_os))

cat("\n─── RECURRENCE-FREE SURVIVAL MODEL ───\n")
print(summary(cox_rfs))

cat("\n─── RESULTS TABLE ───\n")
print(merged_tbl)

cat("\n════════════════════════════════════════════════════════════════════════════\n")
