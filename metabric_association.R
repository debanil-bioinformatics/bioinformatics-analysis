library(tidyverse)
library(stats)

setwd("~/Documents/Cancers")

# ── 1. LOAD DATA ──────────────────────────────────────────────────────────────

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

# Rename columns for easier handling
colnames(clinical_patient) <- gsub(" ", "_", colnames(clinical_patient))
colnames(clinical_sample) <- gsub(" ", "_", colnames(clinical_sample))

clinical <- clinical_sample %>%
  select(`#Patient_Identifier`, 
         Sample_Identifier,
         Tumor_Stage,
         Neoplasm_Histologic_Grade,
         `ER_Status`,
         `PR_Status`,
         `HER2_Status`,
         Tumor_Size) %>%
  rename(PatientID = `#Patient_Identifier`,
         SampleID = Sample_Identifier,
         Stage = Tumor_Stage,
         Grade = Neoplasm_Histologic_Grade,
         ER = `ER_Status`,
         PR = `PR_Status`,
         HER2 = `HER2_Status`) %>%
  left_join(
    clinical_patient %>%
      select(`#Patient_Identifier`,
             `Age_at_Diagnosis`,
             `Lymph_nodes_examined_positive`,
             `Inferred_Menopausal_State`,
             `Pam50_+_Claudin-low_subtype`) %>%
      rename(PatientID = `#Patient_Identifier`,
             Age = `Age_at_Diagnosis`,
             Lymph_Nodes = `Lymph_nodes_examined_positive`,
             Menopausal = `Inferred_Menopausal_State`,
             PAM50 = `Pam50_+_Claudin-low_subtype`),
    by = "PatientID"
  )

# ── 4. CATEGORIZE VARIABLES ───────────────────────────────────────────────────

clinical <- clinical %>%
  mutate(
    Age_Group = case_when(
      Age < 60 ~ "<60",
      Age >= 60 ~ ">=60",
      TRUE ~ NA_character_
    ),
    
    Nodal_Status = case_when(
      Lymph_Nodes == 0 ~ "Negative",
      Lymph_Nodes >= 1 ~ "Positive",
      TRUE ~ NA_character_
    ),
    
    Tumor_Size_Cat = case_when(
      Tumor_Size <= 20 ~ "T1",
      Tumor_Size > 20 & Tumor_Size <= 50 ~ "T2",
      Tumor_Size > 50 ~ "T3",
      TRUE ~ NA_character_
    )
  ) %>%
  select(-Age, -Lymph_Nodes, -Tumor_Size)

# ── 5. MERGE EXPRESSION WITH CLINICAL DATA ────────────────────────────────────

data <- bub1b_expr %>%
  left_join(mad2l1_expr, by = "SampleID") %>%
  left_join(clinical, by = "SampleID") %>%
  drop_na(BUB1B_expr, MAD2L1_expr) %>%
  # Filter out Stage 0 (non-invasive/carcinoma in situ)
  filter(Stage != "0" | is.na(Stage))

cat("Total samples with expression and clinical data:", nrow(data), "\n")

# ── 6. CATEGORIZE BUB1B AND MAD2L1 BY MEDIAN ──────────────────────────────────

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

# ── 7. ASSOCIATION TEST FUNCTION ──────────────────────────────────────────────

compute_association <- function(data, var, display_name = var, gene_status) {
  sub_data <- data %>% 
    filter(!is.na(.data[[var]]), !is.na(.data[[gene_status]]))
  
  # Create contingency table
  table_data <- table(sub_data[[var]], sub_data[[gene_status]])
  table_data <- table_data[rowSums(table_data) > 0, , drop = FALSE]
  
  if (nrow(table_data) == 0) {
    return(list(results = list(), display_name = display_name, p_value = "NA"))
  }
  
  # Calculate percentages
  total_per_level <- rowSums(table_data)
  perc_data <- sweep(table_data, 1, total_per_level, "/") * 100
  
  # Chi-square test
  chi_test <- tryCatch(
    chisq.test(table_data),
    warning = function(w) chisq.test(table_data, simulate.p.value = TRUE)
  )
  
  p_val <- ifelse(chi_test$p.value < 0.0001, 
                  "<0.0001", 
                  sprintf("%.4f", chi_test$p.value))
  
  # Ensure both columns exist
  if (!"Low" %in% colnames(table_data)) {
    table_data <- cbind(table_data, Low = rep(0, nrow(table_data)))
    perc_data <- cbind(perc_data, Low = rep(0, nrow(perc_data)))
  }
  if (!"High" %in% colnames(table_data)) {
    table_data <- cbind(table_data, High = rep(0, nrow(table_data)))
    perc_data <- cbind(perc_data, High = rep(0, nrow(perc_data)))
  }
  
  results <- list()
  for (level in rownames(table_data)) {
    results[[level]] <- c(
      Low = sprintf("%d (%.1f%%)", table_data[level, "Low"], perc_data[level, "Low"]),
      High = sprintf("%d (%.1f%%)", table_data[level, "High"], perc_data[level, "High"]),
      pValue = ""
    )
  }
  
  if (length(results) > 0) {
    results[[names(results)[1]]]["pValue"] <- p_val
  }
  
  return(list(results = results, display_name = display_name, p_value = p_val))
}

# ── 8. CLEAN PAM50 SUBTYPES ──────────────────────────────────────────────────

# Define valid PAM50 subtypes (remove corrupted entries)
valid_pam50 <- c("LumA", "LumB", "Her2", "Basal", "claudin-low", "Normal")

data <- data %>%
  mutate(
    PAM50 = ifelse(PAM50 %in% valid_pam50, PAM50, NA_character_)
  )

# ── 9. VARIABLES TO ANALYZE ──────────────────────────────────────────────────

vars <- list(
  list(var = "Age_Group", display_name = "Age Group"),
  list(var = "Stage", display_name = "Tumor Stage"),
  list(var = "Nodal_Status", display_name = "Nodal Involvement"),
  list(var = "Grade", display_name = "Histologic Grade"),
  list(var = "Menopausal", display_name = "Menopausal State"),
  list(var = "PAM50", display_name = "PAM50 + Claudin-low Subtype")
)

# ── 10. GENERATE ASSOCIATION TABLES ────────────────────────────────────────────

# BUB1B Analysis
results_bub1b <- lapply(vars, function(v) {
  compute_association(data, v$var, v$display_name, "BUB1B_Status")
})

# MAD2L1 Analysis
results_mad2l1 <- lapply(vars, function(v) {
  compute_association(data, v$var, v$display_name, "MAD2L1_Status")
})

# ── 11. FORMAT AND SAVE TABLES AS PDF ────────────────────────────────────────

# Function to create table lines
create_table <- function(results, gene_name, cohort_name = "METABRIC") {
  table_lines <- c(
    sprintf("Supplementary Table: Association between %s expression and clinicopathological features in %s", 
            gene_name, cohort_name),
    "",
    sprintf("%-35s | %-20s | %-20s | %-10s", 
            "Variable", paste0(gene_name, " Low N (%)"), paste0(gene_name, " High N (%)"), "p-Value"),
    paste(rep("-", 90), collapse = "")
  )
  
  for (result in results) {
    # Variable name row
    table_lines <- c(table_lines,
                     sprintf("%-35s | %-20s | %-20s | %-10s",
                             result$display_name, "", "", ""))
    
    # Category rows
    level_names <- names(result$results)
    for (i in seq_along(level_names)) {
      level <- level_names[i]
      line <- sprintf("  %-33s | %-20s | %-20s | %-10s",
                      level,
                      result$results[[level]]["Low"],
                      result$results[[level]]["High"],
                      if (i == 1) result$p_value else "")
      table_lines <- c(table_lines, line)
    }
    table_lines <- c(table_lines, "")
  }
  
  return(table_lines)
}

# Create tables
table_bub1b <- create_table(results_bub1b, "BUB1B")
table_mad2l1 <- create_table(results_mad2l1, "MAD2L1")

# ── 11. SAVE TO PDF ──────────────────────────────────────────────────────────

# BUB1B PDF
pdf("BUB1B_Clinicopath_Association_Metabric_FullCohort.pdf", width = 11, height = 14)
par(mar = c(1, 1, 1, 1))
plot.new()
text(0.05, 0.95, paste(table_bub1b, collapse = "\n"), adj = c(0, 1), family = "mono", cex = 0.65)
dev.off()

# MAD2L1 PDF
pdf("MAD2L1_Clinicopath_Association_Metabric_FullCohort.pdf", width = 11, height = 14)
par(mar = c(1, 1, 1, 1))
plot.new()
text(0.05, 0.95, paste(table_mad2l1, collapse = "\n"), adj = c(0, 1), family = "mono", cex = 0.65)
dev.off()

cat("✅ Analysis complete!\n")
cat("PDFs saved:\n")
cat("  - BUB1B_Clinicopath_Association_Metabric_FullCohort.pdf\n")
cat("  - MAD2L1_Clinicopath_Association_Metabric_FullCohort.pdf\n")
