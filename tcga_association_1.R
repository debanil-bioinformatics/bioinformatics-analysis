library(tidyverse)
library(stats)

# Set your working directory
setwd("~/Documents/Cancers")

cat("TCGA BRCA Association Study: BUB1B and MAD2L1\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

# ── 1. LOAD DATA ──────────────────────────────────────────────────────────────

cat("Loading TCGA BRCA data...\n")
data <- read_csv("TCGA_BRCA.csv", show_col_types = FALSE)

cat("Total samples:", nrow(data), "\n")
cat("Columns:", ncol(data), "\n\n")

# ── 2. EXTRACT AND RENAME KEY COLUMNS ─────────────────────────────────────────

# PAM50 has been removed from this block
data <- data %>%
  select(
    sample,
    BUB1B,
    MAD2L1,
    sample_type,
    Age_at_Initial_Pathologic_Diagnosis_nature2012,
    AJCC_Stage_nature2012,
    Node_Coded_nature2012,
    histological_type,
    menopause_status
  ) %>%
  rename(
    Age = Age_at_Initial_Pathologic_Diagnosis_nature2012,
    Stage_Raw = AJCC_Stage_nature2012,
    Nodal = Node_Coded_nature2012,
    Histology = histological_type,
    Menopause_Raw = menopause_status
  )

# ── 3. PROCESS CLINICAL VARIABLES ─────────────────────────────────────────────

cat("Processing clinical variables...\n\n")

# PAM50 mutation handling has been removed from this block
data <- data %>%
  mutate(
    # Age Group: Binary (<60 vs >=60)
    Age_Group = factor(
      case_when(
        Age < 60  ~ "<60",
        Age >= 60 ~ ">=60",
        TRUE      ~ NA_character_
      ),
      levels = c("<60", ">=60")
    ),
    
    # Tumor Stage: Collapsed strictly into I, II, III, IV
    Stage = factor(
      case_when(
        str_detect(Stage_Raw, "^Stage I[AB]?$") | Stage_Raw == "Stage I"   ~ "I",
        str_detect(Stage_Raw, "^Stage II[AB]?$") | Stage_Raw == "Stage II" ~ "II",
        str_detect(Stage_Raw, "^Stage III[ABC]?$") | Stage_Raw == "Stage III" ~ "III",
        Stage_Raw == "Stage IV"                                            ~ "IV",
        TRUE                                                               ~ NA_character_
      ),
      levels = c("I", "II", "III", "IV")
    ),
    
    # Menopausal State: Shortened strictly to Pre, Peri, Post by prefix matching
    Menopause = factor(
      case_when(
        str_detect(Menopause_Raw, "^Pre")  ~ "Pre",
        str_detect(Menopause_Raw, "^Peri") ~ "Peri",
        str_detect(Menopause_Raw, "^Post") ~ "Post",
        TRUE                               ~ NA_character_ 
      ),
      levels = c("Pre", "Peri", "Post")
    ),
    
    # Keep other molecular and descriptive parameters uncollapsed as requested
    Nodal_Status = factor(Nodal),
    Histology = factor(Histology)
  ) %>%
  select(-Stage_Raw, -Menopause_Raw, -Nodal, -Age)

# Filter for Primary Tumor samples only and clean genomic expression NAs
data_clean <- data %>%
  filter(sample_type == "Primary Tumor") %>%
  drop_na(BUB1B, MAD2L1)

cat("Primary tumor samples:", nrow(data_clean), "\n\n")

# ── 4. CATEGORIZE BUB1B AND MAD2L1 BY MEDIAN ──────────────────────────────────

data_clean <- data_clean %>%
  mutate(
    BUB1B_Status = factor(
      ifelse(BUB1B >= median(BUB1B, na.rm = TRUE), "High", "Low"),
      levels = c("Low", "High")
    ),
    MAD2L1_Status = factor(
      ifelse(MAD2L1 >= median(MAD2L1, na.rm = TRUE), "High", "Low"),
      levels = c("Low", "High")
    )
  )

cat("BUB1B median:", median(data_clean$BUB1B, na.rm = TRUE), "\n")
cat("MAD2L1 median:", median(data_clean$MAD2L1, na.rm = TRUE), "\n\n")

# ── 5. ASSOCIATION TEST FUNCTION ──────────────────────────────────────────────

compute_association <- function(data, var, display_name = var, gene_status) {
  sub_data <- data %>% 
    filter(!is.na(.data[[var]]), !is.na(.data[[gene_status]]))
  
  table_data <- table(sub_data[[var]], sub_data[[gene_status]])
  table_data <- table_data[rowSums(table_data) > 0, , drop = FALSE]
  
  if (nrow(table_data) == 0) {
    return(list(results = list(), display_name = display_name, p_value = "NA"))
  }
  
  total_per_level <- rowSums(table_data)
  perc_data <- sweep(table_data, 1, total_per_level, "/") * 100
  
  chi_test <- tryCatch(
    chisq.test(table_data),
    warning = function(w) chisq.test(table_data, simulate.p.value = TRUE)
  )
  
  p_val <- ifelse(chi_test$p.value < 0.0001, 
                  "<0.0001", 
                  sprintf("%.4f", chi_test$p.value))
  
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

# ── 6. VARIABLES TO ANALYZE ──────────────────────────────────────────────────

# PAM50 Subtype row has been removed from this list
vars <- list(
  list(var = "Age_Group", display_name = "Age Group"),
  list(var = "Stage", display_name = "Tumor Stage"),
  list(var = "Nodal_Status", display_name = "Nodal Involvement"),
  list(var = "Histology", display_name = "Histological Type"),
  list(var = "Menopause", display_name = "Menopausal State")
)

# ── 7. RUN ASSOCIATION TESTS ──────────────────────────────────────────────────

cat("Running chi-square association tests...\n\n")

results_bub1b <- lapply(vars, function(v) {
  compute_association(data_clean, v$var, v$display_name, "BUB1B_Status")
})

results_mad2l1 <- lapply(vars, function(v) {
  compute_association(data_clean, v$var, v$display_name, "MAD2L1_Status")
})

# ── 8. FORMAT AND SAVE TABLES AS PDF ──────────────────────────────────────────

create_table <- function(results, gene_name, cohort_name = "TCGA BRCA") {
  table_lines <- c(
    sprintf("Supplementary Table: Association between %s expression and clinicopathological features in %s", 
            gene_name, cohort_name),
    "",
    sprintf("%-45s | %-20s | %-20s | %-10s", 
            "Variable", paste0(gene_name, " Low N (%)"), paste0(gene_name, " High N (%)"), "p-Value"),
    paste(rep("-", 105), collapse = "")
  )
  
  for (result in results) {
    table_lines <- c(table_lines,
                     sprintf("%-45s | %-20s | %-20s | %-10s",
                             result$display_name, "", "", ""))
    
    level_names <- names(result$results)
    for (i in seq_along(level_names)) {
      level <- level_names[i]
      level_trimmed <- substr(level, 1, 40)
      line <- sprintf("  %-43s | %-20s | %-20s | %-10s",
                      level_trimmed,
                      result$results[[level]]["Low"],
                      result$results[[level]]["High"],
                      if (i == 1) result$p_value else "")
      table_lines <- c(table_lines, line)
    }
    table_lines <- c(table_lines, "")
  }
  
  return(table_lines)
}

table_bub1b <- create_table(results_bub1b, "BUB1B")
table_mad2l1 <- create_table(results_mad2l1, "MAD2L1")

# ── 9. SAVE TO PDF ──────────────────────────────────────────────────────────

cat("Creating PDF tables...\n\n")

pdf("BUB1B_Clinicopath_Association_TCGA_FullCohort.pdf", width = 12, height = 11)
par(mar = c(1, 1, 1, 1))
plot.new()
text(0.02, 0.98, paste(table_bub1b, collapse = "\n"), adj = c(0, 1), family = "mono", cex = 0.6)
dev.off()

pdf("MAD2L1_Clinicopath_Association_TCGA_FullCohort.pdf", width = 12, height = 11)
par(mar = c(1, 1, 1, 1))
plot.new()
text(0.02, 0.98, paste(table_mad2l1, collapse = "\n"), adj = c(0, 1), family = "mono", cex = 0.6)
dev.off()

cat("✅ TCGA Association Study Complete (Excluding PAM50)!\n")