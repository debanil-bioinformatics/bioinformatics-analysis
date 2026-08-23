library(survival)
library(broom)
library(dplyr)
library(tidyr)
library(flextable)
library(officer)
library(readr)

# Set working directory
setwd("~/Documents/Cancers")

cat("Starting TCGA BRCA Multivariate Cox Regression Analysis...\n\n")

# ── 1. LOAD DATA ──────────────────────────────────────────────────────────────

cat("Loading TCGA BRCA data...\n")
data <- read_csv("TCGA_BRCA.csv", show_col_types = FALSE)

cat("Total samples:", nrow(data), "\n")
cat("Columns:", ncol(data), "\n\n")

# ── 2. PROCESS AND BINARIZE CLINICAL COVARIATES ────────────────────────────────

cat("Processing and binarizing variables...\n\n")

data_processed <- data %>%
  filter(sample_type == "Primary Tumor") %>%
  mutate(
    # Age Group: <60 vs >=60
    Age_Group = factor(
      case_when(
        Age_at_Initial_Pathologic_Diagnosis_nature2012 < 60  ~ "Low (<60)",
        Age_at_Initial_Pathologic_Diagnosis_nature2012 >= 60 ~ "High (>=60)",
        TRUE ~ NA_character_
      ),
      levels = c("Low (<60)", "High (>=60)")
    ),
    
    # Tumor Stage: Low (I-II) vs High (III-IV)
    Stage = factor(
      case_when(
        str_detect(AJCC_Stage_nature2012, "^Stage I[AB]?$") | AJCC_Stage_nature2012 == "Stage I"   ~ "Low (I-II)",
        str_detect(AJCC_Stage_nature2012, "^Stage II[AB]?$") | AJCC_Stage_nature2012 == "Stage II" ~ "Low (I-II)",
        str_detect(AJCC_Stage_nature2012, "^Stage III[ABC]?$") | AJCC_Stage_nature2012 == "Stage III" ~ "High (III-IV)",
        AJCC_Stage_nature2012 == "Stage IV"                                                    ~ "High (III-IV)",
        TRUE ~ NA_character_
      ),
      levels = c("Low (I-II)", "High (III-IV)")
    ),
    
    # Nodal Status: Negative vs Positive
    Nodal_Status = factor(
      case_when(
        Node_Coded_nature2012 == "Negative" ~ "Negative",
        Node_Coded_nature2012 == "Positive" ~ "Positive",
        TRUE ~ NA_character_
      ),
      levels = c("Negative", "Positive")
    ),
    
    # Histological Type: Infiltrating Ductal vs Others
    Histology = factor(
      case_when(
        histological_type == "Infiltrating Ductal Carcinoma" ~ "Infiltrating Ductal",
        !is.na(histological_type)                             ~ "Other Histologies",
        TRUE ~ NA_character_
      ),
      levels = c("Infiltrating Ductal", "Other Histologies")
    ),
    
    # Menopausal State: Pre/Peri vs Post
    Menopause = factor(
      case_when(
        str_detect(menopause_status, "^Pre") | str_detect(menopause_status, "^Peri") ~ "Pre/Peri",
        str_detect(menopause_status, "^Post")                                       ~ "Post",
        TRUE ~ NA_character_
      ),
      levels = c("Pre/Peri", "Post")
    ),
    
    # PAM50 Subtype: Luminal vs Non-Luminal
    PAM50_Group = factor(
      case_when(
        PAM50_mRNA_nature2012 %in% c("Luminal A", "Luminal B") ~ "Luminal",
        PAM50_mRNA_nature2012 %in% c("Basal-like", "HER2-enriched", "Normal-like") ~ "Non-Luminal",
        TRUE ~ NA_character_
      ),
      levels = c("Luminal", "Non-Luminal")
    ),
    
    # Target Genes: Binarized by primary tumor cohort median expression values
    BUB1B_Status = factor(
      ifelse(BUB1B >= median(BUB1B, na.rm = TRUE), "High", "Low"),
      levels = c("Low", "High")
    ),
    MAD2L1_Status = factor(
      ifelse(MAD2L1 >= median(MAD2L1, na.rm = TRUE), "High", "Low"),
      levels = c("Low", "High")
    )
  ) %>%
  # Filter missing baseline elements for survival endpoints
  drop_na(OS_event_nature2012, OS_Time_nature2012, DFI, DFI.time, BUB1B, MAD2L1)

# ── 3. DEFINE SURVIVAL OUTCOMES AND RUN MODELS ────────────────────────────────

# Formulate formulas exactly as requested using both target markers jointly
os_formula  <- Surv(OS_Time_nature2012, OS_event_nature2012) ~ Age_Group + Stage + Nodal_Status + Histology + Menopause + PAM50_Group + BUB1B_Status + MAD2L1_Status
dfi_formula <- Surv(DFI.time, DFI) ~ Age_Group + Stage + Nodal_Status + Histology + Menopause + PAM50_Group + BUB1B_Status + MAD2L1_Status

cat("Fitting Multivariate Cox proportional hazards models...\n")
cox_os  <- coxph(os_formula, data = data_processed)
cox_dfi <- coxph(dfi_formula, data = data_processed)

# ── 4. EXTRACT HARVESTED METRICS ──────────────────────────────────────────────

clean_cox_results <- function(model_fit) {
  tidy(model_fit, exponentiate = TRUE, conf.int = TRUE) %>%
    select(term, estimate, conf.low, conf.high, p.value)
}

res_os  <- clean_cox_results(cox_os)
res_dfi <- clean_cox_results(cox_dfi)

# ── 5. MERGE ENDPOINTS & BUILD STRUCTURE ──────────────────────────────────────

# Set up visual naming mapping
display_mapping <- data.frame(
  term = c("Age_GroupHigh (>=60)", "StageHigh (III-IV)", "Nodal_StatusPositive", 
           "HistologyOther Histologies", "MenopausePost", "PAM50_GroupNon-Liminal", 
           "BUB1B_StatusHigh", "MAD2L1_StatusHigh"),
  Variable = c("Age Group", "Tumor Stage", "Nodal Involvement", 
               "Histological Type", "Menopausal State", "PAM50 Subtype", 
               "BUB1B Status", "MAD2L1 Status"),
  Category = c("High (>=60)", "High (III-IV)", "Positive", 
               "Other Histologies", "Post", "Non-Luminal", 
               "High", "High"),
  stringsAsFactors = FALSE
)

# Reference definitions matching your preferred publication style
reference_rows <- data.frame(
  Variable = c("Age Group", "Tumor Stage", "Nodal Involvement", "Histological Type", "Menopausal State", "PAM50 Subtype", "BUB1B Status", "MAD2L1 Status"),
  Category = c("Low (<60)", "Low (I-II)", "Negative", "Infiltrating Ductal", "Pre/Peri", "Luminal", "Low", "Low"),
  OS_HR = "1.00", OS_P = "", DFI_HR = "1.00", DFI_P = "",
  stringsAsFactors = FALSE
)

final_data <- display_mapping %>%
  left_join(res_os %>% select(term, est_os = estimate, low_os = conf.low, high_os = conf.high, p_os = p.value), by = "term") %>%
  left_join(res_dfi %>% select(term, est_dfi = estimate, low_dfi = conf.low, high_dfi = conf.high, p_dfi = p.value), by = "term") %>%
  mutate(
    OS_HR = sprintf("%.2f (%.2f-%.2f)", est_os, low_os, high_os),
    OS_P  = ifelse(p_os < 0.001, "<0.001", sprintf("%.3f", p_os)),
    DFI_HR = sprintf("%.2f (%.2f-%.2f)", est_dfi, low_dfi, high_dfi),
    DFI_P  = ifelse(p_dfi < 0.001, "<0.001", sprintf("%.3f", p_dfi))
  ) %>%
  select(Variable, Category, OS_HR, OS_P, DFI_HR, DFI_P)

# Combine references and outcomes hierarchically
structured_table <- data.frame()
for(i in 1:nrow(reference_rows)) {
  structured_table <- rbind(structured_table, reference_rows[i, ])
  match_row <- final_data %>% filter(Variable == reference_rows$Variable[i])
  if(nrow(match_row) > 0) {
    # Clean redundant parent variable text to achieve a minimal presentation look
    match_row$Variable <- "" 
    structured_table <- rbind(structured_table, match_row)
  }
}

# ── 6. GENERATE PUBLICATION-READY FLEXTABLE ───────────────────────────────────

cat("Building output table layouts...\n")

ft <- flextable(structured_table) %>%
  set_header_labels(
    Variable = "Variables", Category = "Category",
    OS_HR = "Overall Survival\nHazard Ratio (95% CI)", OS_P = "Overall Survival\np-Value",
    DFI_HR = "DFI\nHazard Ratio (95% CI)", DFI_P = "DFI\np-Value"
  ) %>%
  align(align = "center", part = "header") %>%
  align(j = 3:6, align = "center", part = "body") %>%
  font(fontname = "Arial", part = "all") %>%
  fontsize(size = 10, part = "all") %>%
  border_remove() %>%
  hline_top(part = "header", border = fp_border(color = "black", width = 2)) %>%
  hline_bottom(part = "body", border = fp_border(color = "black", width = 2)) %>%
  autofit() %>%
  fit_to_width(11)

# ── 7. SAVE OUTPUT EXPORTS ───────────────────────────────────────────────────

cat("Saving multivariate Cox regression artifacts to disk...\n")

# Save high-resolution imagery and compiled PDF table
save_as_image(ft, path = "tcga_cox_results.png", zoom = 2)

pdf("TCGA_Multivariate_Cox_Regression_BUB1B_MAD2L1.pdf", width = 12, height = 8)
grid.newpage()
img <- readPNG("tcga_cox_results.png")
grid.raster(img)
dev.off()

cat("\n════════════════════════════════════════════════════════════════════════════\n")
cat("✅ MULTIVARIATE COX REGRESSION FOR TCGA COMPLETE!\n")
cat("Artifact generated: TCGA_Multivariate_Cox_Regression_BUB1B_MAD2L1.pdf\n")
cat("════════════════════════════════════════════════════════════════════════════\n")