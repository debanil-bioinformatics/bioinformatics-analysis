############################################################
# Author: Debanil Dhar
# Year: 2025
#
# Description:
# Generalized clinicopathological association analysis
# across TCGA projects. Computes associations between
# a binary molecular status variable and clinical
# features using chi-square tests and generates a
# publication-ready supplementary table (PDF).
#
# Input:
# - TCGA project–specific CSV with clinical variables
#   and a binary grouping column (e.g., Status)
#
# Output:
# - PDF table summarizing associations and p-values
#
# Notes:
# - Designed for reuse across TCGA cancer types
# - Gracefully handles missing clinical variables
############################################################

## Load required libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stats)
})

## ---------------- USER PARAMETERS ----------------
tcga_project  <- "TCGA-BRCA"
input_file   <- "tcga_clinical_data.csv"
status_col   <- "Status"          # Binary grouping variable
group_levels <- c("Both.High", "Others")
output_dir   <- "results"
## -------------------------------------------------

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

## Load data
data <- read.csv(input_file, na.strings = c("", "NA"), stringsAsFactors = FALSE)

stopifnot(status_col %in% colnames(data))

## Ensure grouping variable is factored
data[[status_col]] <- factor(data[[status_col]], levels = group_levels)

## Recode clinical variables (only if present)
data <- data %>%
  mutate(
    Age_Group = if ("Age.at.Diagnosis" %in% colnames(.))
      factor(ifelse(Age.at.Diagnosis < 60, "<60", ">=60"),
             levels = c("<60", ">=60")) else NA,
    
    Tumor_Stage = if ("clinical_stage" %in% colnames(.))
      case_when(
        clinical_stage %in% c("Stage I", "Stage IA", "Stage IB") ~ "Stage I",
        clinical_stage %in% c("Stage II", "Stage IIA", "Stage IIB") ~ "Stage II",
        clinical_stage %in% c("Stage III", "Stage IIIA", "Stage IIIB", "Stage IIIC") ~ "Stage III",
        clinical_stage %in% c("Stage IV") ~ "Stage IV",
        TRUE ~ NA_character_
      ) else NA,
    
    Nodal_Status = if ("clinical_N" %in% colnames(.))
      case_when(
        clinical_N == "N0" ~ "N0",
        clinical_N %in% c("N1", "N1a", "N1b", "N1c") ~ "N1",
        clinical_N %in% c("N2", "N2a") ~ "N2",
        clinical_N %in% c("N3", "N3a", "N3b", "N3c") ~ "N3",
        TRUE ~ NA_character_
      ) else NA,
    
    Metastatic_Status = if ("clinical_M" %in% colnames(.))
      case_when(
        clinical_M %in% c("M0", "cM0") ~ "M0",
        clinical_M == "M1" ~ "M1",
        TRUE ~ NA_character_
      ) else NA,
    
    Tumor_Size = if ("pathologic_T" %in% colnames(.))
      case_when(
        pathologic_T %in% c("T1", "T1a", "T1b", "T1c") ~ "T1",
        pathologic_T %in% c("T2", "T2b") ~ "T2",
        pathologic_T == "T3" ~ "T3",
        pathologic_T %in% c("T4b", "T4d") ~ "T4",
        TRUE ~ NA_character_
      ) else NA
  )

## Association test function
compute_association <- function(data, var, display_name) {
  
  if (!var %in% colnames(data)) {
    return(list(results = list(), display_name = display_name, p_value = NA))
  }
  
  sub_data <- data %>% filter(!is.na(.data[[var]]), !is.na(.data[[status_col]]))
  table_data <- table(sub_data[[var]], sub_data[[status_col]])
  
  if (nrow(table_data) == 0 || ncol(table_data) < 2) {
    return(list(results = list(), display_name = display_name, p_value = NA))
  }
  
  chi_test <- suppressWarnings(
    chisq.test(table_data, simulate.p.value = TRUE)
  )
  
  p_val <- ifelse(chi_test$p.value < 0.0001,
                  "<0.0001",
                  sprintf("%.4f", chi_test$p.value))
  
  perc_data <- sweep(table_data, 1, rowSums(table_data), "/") * 100
  
  results <- lapply(rownames(table_data), function(level) {
    c(
      table_data[level, group_levels[2]],
      sprintf("%.1f", perc_data[level, group_levels[2]]),
      table_data[level, group_levels[1]],
      sprintf("%.1f", perc_data[level, group_levels[1]])
    )
  })
  
  names(results) <- rownames(table_data)
  
  list(results = results, display_name = display_name, p_value = p_val)
}

## Variables to test
variables <- list(
  c("Age_Group", "Age Group"),
  c("Tumor_Stage", "Tumor Stage"),
  c("Nodal_Status", "Nodal Involvement"),
  c("Metastatic_Status", "Metastatic Status"),
  c("Tumor_Size", "Tumor Size")
)

results <- lapply(variables, function(v)
  compute_association(data, v[1], v[2]))

## Generate PDF
pdf_file <- file.path(
  output_dir,
  paste0(tcga_project, "_clinicopathological_association.pdf")
)

pdf(pdf_file, width = 10, height = 14)
par(mar = c(2, 2, 2, 2))
plot.new()

text(
  0.05, 0.95,
  paste(
    c(
      paste("Supplementary Table:",
            "Association between molecular status and clinicopathological features in",
            tcga_project),
      "",
      unlist(lapply(results, function(r)
        c(r$display_name, paste("p-value:", r$p_value), "")))
    ),
    collapse = "\n"
  ),
  adj = c(0, 1),
  family = "mono",
  cex = 0.7
)

dev.off()

message("Clinicopathological association analysis completed for ", tcga_project)
