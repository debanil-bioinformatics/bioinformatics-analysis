library(dplyr)
library(stats)

setwd("~/Documents/TCGA Analysis")

## --------------------------------------------------
## 1. Read cleaned dataset
## --------------------------------------------------

data <- read.csv(
  "TCGA HR NEG.csv",
  stringsAsFactors = FALSE,
  na.strings = c("", "NA")
)

cat("Initial sample size (all samples):", nrow(data), "\n")

## --------------------------------------------------
## 1.1 Tumour-only filtering
## --------------------------------------------------

if (!"sample_type" %in% colnames(data)) {
  stop("ERROR: 'sample_type' column not found. Cannot filter tumour samples.")
}

data <- data %>%
  filter(sample_type == "Tumour")

cat("Final sample size (tumour only):", nrow(data), "\n")

if (nrow(data) == 0) {
  stop("ERROR: Tumour filtering removed all samples.")
}

## --------------------------------------------------
## 2. Define derived clinical variables
## --------------------------------------------------

data <- data %>%
  mutate(
    
    ## Age group
    Age_Group = factor(
      ifelse(age_at_initial_pathologic_diagnosis < 40, "<40", ">=40"),
      levels = c("<40", ">=40")
    ),
    
    ## Histology
    Histology = factor(
      ifelse(grepl("^8500", icd_o_3_histology), "IDC", "Others"),
      levels = c("IDC", "Others")
    ),
    
    ## Tumour size (T): T1–T4
    Tumor_size = factor(
      case_when(
        grepl("^T1", Tumor_nature2012) ~ "T1",
        grepl("^T2", Tumor_nature2012) ~ "T2",
        grepl("^T3", Tumor_nature2012) ~ "T3",
        grepl("^T4", Tumor_nature2012) ~ "T4",
        TRUE ~ NA_character_
      ),
      levels = c("T1", "T2", "T3", "T4")
    ),
    
    ## Nodal status (N)
    Nodal_Status = factor(
      case_when(
        grepl("^N0", Node_nature2012) ~ "N0",
        grepl("^N1", Node_nature2012) ~ "N1",
        grepl("^N2", Node_nature2012) ~ "N2",
        grepl("^N3", Node_nature2012) ~ "N3",
        TRUE ~ NA_character_
      ),
      levels = c("N0", "N1", "N2", "N3")
    ),
    
    ## Pathologic stage (I–IV; substages collapsed)
    Stage = factor(
      case_when(
        grepl("^Stage IV",  AJCC_Stage_nature2012) ~ "IV",
        grepl("^Stage III", AJCC_Stage_nature2012) ~ "III",
        grepl("^Stage II",  AJCC_Stage_nature2012) ~ "II",
        grepl("^Stage I",   AJCC_Stage_nature2012) ~ "I",
        TRUE ~ NA_character_
      ),
      levels = c("I", "II", "III", "IV")
    ),
    
    ## Menopausal status
    Menopause = factor(
      case_when(
        grepl("^Pre",  menopause_status)  ~ "Pre",
        grepl("^Peri", menopause_status)  ~ "Peri",
        grepl("^Post", menopause_status)  ~ "Post",
        TRUE                              ~ NA_character_
      ),
      levels = c("Pre", "Peri", "Post")
    ),
    
    ## Molecular subtype (curated upstream)
    Subtype = factor(
      HR_status_strat,
      levels = c("HR+", "HER2+", "TNBC", "HR-")
    )
  )

## --------------------------------------------------
## 3. Define analysis-eligible cohort (OPTION 1)
## --------------------------------------------------

analysis_vars <- c(
  "Age_Group",
  "Subtype",
  "Histology",
  "Nodal_Status",
  "Tumor_size",
  "Stage",
  "Menopause"
)

eligible_for_analysis <- data %>%
  filter(
    !is.na(CUEDC2) &
      rowSums(!is.na(select(., all_of(analysis_vars)))) > 0
  )

cat(
  "Samples eligible for median calculation:",
  nrow(eligible_for_analysis), "\n"
)

## --------------------------------------------------
## 4. Median splits (ELIGIBLE SAMPLES ONLY)
## --------------------------------------------------

cuedc2_median <- median(eligible_for_analysis$CUEDC2, na.rm = TRUE)
MKI67_median  <- median(eligible_for_analysis$MKI67,  na.rm = TRUE)

data <- data %>%
  mutate(
    Status = factor(
      ifelse(CUEDC2 >= cuedc2_median, "High", "Low"),
      levels = c("Low", "High")
    ),
    
    MKI67_Group = factor(
      ifelse(MKI67 >= MKI67_median, "High", "Low"),
      levels = c("Low", "High")
    )
  )

## --------------------------------------------------
## 5. Association test function
##    (Chi-square or Fisher’s exact test as appropriate)
## --------------------------------------------------

compute_association <- function(data, var, display_name = var) {
  
  valid <- !is.na(data[[var]]) & !is.na(data$Status)
  
  tab <- table(data[[var]][valid], data$Status[valid])
  tab <- tab[rowSums(tab) > 0, , drop = FALSE]
  
  expected <- suppressWarnings(chisq.test(tab)$expected)
  use_fisher <- any(expected < 5)
  
  test <- if (use_fisher) fisher.test(tab) else chisq.test(tab)
  
  p_val <- ifelse(
    test$p.value < 0.0001,
    "<0.0001",
    sprintf("%.4f", test$p.value)
  )
  
  perc <- sweep(tab, 1, rowSums(tab), "/") * 100
  
  res <- list()
  for (lvl in rownames(tab)) {
    res[[lvl]] <- c(
      Low  = sprintf("%d (%.1f)", tab[lvl, "Low"],  perc[lvl, "Low"]),
      High = sprintf("%d (%.1f)", tab[lvl, "High"], perc[lvl, "High"]),
      pValue = ""
    )
  }
  
  res[[1]]["pValue"] <- p_val
  
  list(display_name = display_name, results = res)
}

## --------------------------------------------------
## 6. Variables for association
## --------------------------------------------------

vars <- list(
  list("Age_Group",    "Age"),
  list("Subtype",      "Molecular subtype"),
  list("MKI67_Group",  "MKI67"),
  list("Histology",    "Histology"),
  list("Nodal_Status", "Nodal status"),
  list("Tumor_size",   "Tumour size (T)"),
  list("Stage",        "Pathologic stage"),
  list("Menopause",    "Menopausal status")
)

results <- lapply(vars, function(v)
  compute_association(data, v[[1]], v[[2]])
)

## --------------------------------------------------
## 7. Build publication-ready PDF table
## --------------------------------------------------

table_lines <- c(
  "Supplementary Table: Association between CUEDC2 expression and clinicopathological features in TCGA BRCA (tumour samples only)",
  "",
  sprintf("%-30s | %-20s | %-20s | %-10s",
          "Variables", "Low CUEDC2 N (%)", "High CUEDC2 N (%)", "p-Value"),
  paste(rep("-", 90), collapse = "")
)

for (r in results) {
  
  table_lines <- c(
    table_lines,
    sprintf("%-30s | %-20s | %-20s | %-10s",
            r$display_name, "", "", "")
  )
  
  for (lvl in names(r$results)) {
    x <- r$results[[lvl]]
    table_lines <- c(
      table_lines,
      sprintf("  %-28s | %-20s | %-20s | %-10s",
              lvl, x["Low"], x["High"], x["pValue"])
    )
  }
  
  table_lines <- c(table_lines, "")
}

encoding_note <- c(
  "",
  "Encoding legend:",
  "- Dataset includes TCGA-BRCA tumour samples only",
  "- CUEDC2 and MKI67 dichotomised using medians calculated from analysis-eligible tumour samples",
  "- Associations tested using Chi-square or Fisher’s exact test as appropriate",
  "- Molecular subtype based on curated HR_status_strat: HR+, HER2+, TNBC, HR-",
  "- Age: <40 vs >=40",
  "- Histology: IDC vs Others",
  "- Tumour size: T1, T2, T3, T4",
  "- Nodal status: N0, N1, N2, N3",
  "- AJCC_Stage_nature2012: I (I/IA/IB), II (II/IIA/IIB), III (III/IIIA/IIIB/IIIC), IV",
  "- Menopause: Pre, Peri, Post"
)

pdf("updated_association_full_cohort.pdf", width = 10, height = 14)
par(mar = c(2, 2, 2, 2))
plot.new()
text(
  0.05, 0.95,
  paste(c(table_lines, encoding_note), collapse = "\n"),
  adj = c(0, 1),
  family = "mono",
  cex = 0.7
)
dev.off()

cat(
  "Association analysis completed successfully ",
  "(tumour-only; medians based on analysis-eligible samples).\n"
)
