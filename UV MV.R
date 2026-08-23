## --------------------------------------------------
## HR− cohort Cox regression analysis (OS & PFI)
## Univariate + Multivariate models
## Factors derived from CODE 1
## CUEDC2 stratified within HR− cohort (median)
## --------------------------------------------------

setwd("~/Documents/TCGA Analysis")

library(dplyr)
library(survival)
library(tidyr)

## --------------------------------------------------
## A. Read annotated dataset
## --------------------------------------------------

tcga <- read.csv("TCGA HR NEG.csv", stringsAsFactors = FALSE)

## HR− tumour-only cohort
tcga_hrneg <- tcga %>%
  filter(sample_type == "Tumour",
         HR_status == "HR-",
         !is.na(CUEDC2))

## --------------------------------------------------
## B. Recreate CODE 1 clinical factors
## --------------------------------------------------

tcga_hrneg <- tcga_hrneg %>%
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
    
    ## Tumour size (T)
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
    
    ## Nodal status
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
    
    ## Stage
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
    
    ## Menopause
    Menopause = factor(
      case_when(
        grepl("^Pre",  menopause_status)  ~ "Pre",
        grepl("^Peri", menopause_status)  ~ "Peri",
        grepl("^Post", menopause_status)  ~ "Post",
        TRUE ~ NA_character_
      ),
      levels = c("Pre", "Peri", "Post")
    )
  )

## --------------------------------------------------
## C. CUEDC2 median split (HR− only)
## --------------------------------------------------

median_cuedc2_hrneg <- median(tcga_hrneg$CUEDC2, na.rm = TRUE)

tcga_hrneg <- tcga_hrneg %>%
  mutate(
    CUEDC2_group = factor(ifelse(CUEDC2 >= median_cuedc2_hrneg, "High", "Low"),
                          levels = c("Low", "High"))
  )

## --------------------------------------------------
## D. Time variables (MONTHS)
## --------------------------------------------------

tcga_hrneg <- tcga_hrneg %>%
  mutate(
    OS_months  = OS.time  / 30.44,
    PFI_months = PFI.time / 30.44
  )

## --------------------------------------------------
## E. Variables for models (from CODE 1)
## --------------------------------------------------

cox_vars <- c(
  "CUEDC2_group",
  "Age_Group",
  "Histology",
  "Tumor_size",
  "Nodal_Status",
  "Stage",
  "Menopause"
)

## --------------------------------------------------
## F. UNIVARIATE COX
## --------------------------------------------------

run_univariate <- function(data, time, event, vars) {
  res <- lapply(vars, function(v) {
    df <- data %>% select(all_of(c(time, event, v))) %>% drop_na()
    if (nrow(df) < 20) return(NULL)
    
    f <- as.formula(paste0("Surv(", time, ",", event, ") ~ ", v))
    fit <- coxph(f, data = df)
    s <- summary(fit)
    
    data.frame(
      Variable = v,
      Level = rownames(s$coefficients),
      HR = round(s$coefficients[,"exp(coef)"], 3),
      CI_low = round(s$conf.int[,"lower .95"], 3),
      CI_high = round(s$conf.int[,"upper .95"], 3),
      P_value = signif(s$coefficients[,"Pr(>|z|)"], 3),
      stringsAsFactors = FALSE
    )
  })
  
  bind_rows(res)
}

## OS univariate
uni_os <- run_univariate(tcga_hrneg %>% filter(!is.na(OS), !is.na(OS_months)),
                         "OS_months", "OS", cox_vars)

## PFI univariate
uni_pfi <- run_univariate(tcga_hrneg %>% filter(!is.na(PFI), !is.na(PFI_months)),
                          "PFI_months", "PFI", cox_vars)

## --------------------------------------------------
## G. MULTIVARIATE COX
## --------------------------------------------------

## OS multivariate
multi_os_data <- tcga_hrneg %>%
  select(OS_months, OS, all_of(cox_vars)) %>%
  drop_na()

cox_multi_os <- coxph(
  Surv(OS_months, OS) ~ CUEDC2_group + Age_Group + Histology + Tumor_size +
    Nodal_Status + Stage + Menopause,
  data = multi_os_data
)

sum_multi_os <- summary(cox_multi_os)

multi_os <- data.frame(
  Variable = rownames(sum_multi_os$coefficients),
  HR = round(sum_multi_os$coefficients[,"exp(coef)"], 3),
  CI_low = round(sum_multi_os$conf.int[,"lower .95"], 3),
  CI_high = round(sum_multi_os$conf.int[,"upper .95"], 3),
  P_value = signif(sum_multi_os$coefficients[,"Pr(>|z|)"], 3),
  stringsAsFactors = FALSE
)

## PFI multivariate
multi_pfi_data <- tcga_hrneg %>%
  select(PFI_months, PFI, all_of(cox_vars)) %>%
  drop_na()

cox_multi_pfi <- coxph(
  Surv(PFI_months, PFI) ~ CUEDC2_group + Age_Group + Histology + Tumor_size +
    Nodal_Status + Stage + Menopause,
  data = multi_pfi_data
)

sum_multi_pfi <- summary(cox_multi_pfi)

multi_pfi <- data.frame(
  Variable = rownames(sum_multi_pfi$coefficients),
  HR = round(sum_multi_pfi$coefficients[,"exp(coef)"], 3),
  CI_low = round(sum_multi_pfi$conf.int[,"lower .95"], 3),
  CI_high = round(sum_multi_pfi$conf.int[,"upper .95"], 3),
  P_value = signif(sum_multi_pfi$coefficients[,"Pr(>|z|)"], 3),
  stringsAsFactors = FALSE
)

## --------------------------------------------------
## H. EXPORT RESULTS
## --------------------------------------------------

write.csv(uni_os,   "HRneg_CUEDC2_Univariate_OS.csv", row.names = FALSE)
write.csv(uni_pfi,  "HRneg_CUEDC2_Univariate_PFI.csv", row.names = FALSE)
write.csv(multi_os, "HRneg_CUEDC2_Multivariate_OS.csv", row.names = FALSE)
write.csv(multi_pfi,"HRneg_CUEDC2_Multivariate_PFI.csv", row.names = FALSE)

cat("\nCox regression completed:\n",
    "- Univariate OS & PFI\n",
    "- Multivariate OS & PFI\n",
    "Files exported successfully.\n")
