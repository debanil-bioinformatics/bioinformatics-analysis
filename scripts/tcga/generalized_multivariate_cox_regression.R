############################################################
# Author: Debanil Dhar
# Year: 2025
#
# Description:
# Generalized multivariate Cox proportional hazards
# regression framework for survival analysis across
# cancer cohorts (TCGA, METABRIC, or custom datasets).
#
# Input:
# - CSV file with survival outcomes, a binary biomarker,
#   and clinical covariates
#
# Output:
# - Publication-ready PDF with Cox regression tables
#
# Notes:
# - Supports multiple survival endpoints (e.g., OS, RFS)
# - Designed for reuse across datasets and biomarkers
############################################################

## Load required libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(survival)
  library(gtsummary)
  library(grid)
  library(officer)
  library(flextable)
  library(png)
})

## ---------------- USER PARAMETERS ----------------

# Input
input_file <- "clinical_survival_data.csv"

# Biomarker / grouping variable
status_col    <- "Status"          # e.g., High / Low
status_levels <- c("Low", "High")  # reference first

# Survival endpoints (define as list)
survival_endpoints <- list(
  OS  = list(time = "OS_time",  event = "OS_event"),
  RFS = list(time = "RFS_time", event = "RFS_event")
)

# Covariates to include in the model
covariates <- c(
  "Age_Group",
  "Tumor_Stage_Binary",
  "Nodal_Status",
  "Tumor_Size_Binary",
  "Grade_Binary",
  "Menopausal_Binary"
)

# Output
output_pdf <- "multivariate_cox_regression_results.pdf"

## -------------------------------------------------

## Load data
data <- read.csv(input_file, na.strings = c("", "NA"), stringsAsFactors = FALSE)

## Sanity checks
stopifnot(status_col %in% colnames(data))
stopifnot(all(unlist(lapply(survival_endpoints, unlist)) %in% colnames(data)))

## Factorize biomarker
data[[status_col]] <- factor(data[[status_col]], levels = status_levels)

## Select complete cases
model_vars <- c(
  status_col,
  covariates,
  unlist(lapply(survival_endpoints, unlist))
)

complete_data <- data %>%
  select(all_of(model_vars)) %>%
  na.omit()

message("Total samples: ", nrow(data))
message("Samples with complete data: ", nrow(complete_data))

## Fit Cox models
cox_tables <- lapply(names(survival_endpoints), function(endpoint) {
  
  ep <- survival_endpoints[[endpoint]]
  
  formula_str <- paste(
    "Surv(", ep$time, ",", ep$event, ") ~",
    paste(c(status_col, covariates), collapse = " + ")
  )
  
  cox_model <- coxph(as.formula(formula_str), data = complete_data)
  
  tbl_regression(cox_model, exponentiate = TRUE) %>%
    modify_header(label = "**Variable**") %>%
    modify_caption(paste(endpoint, "Cox Model"))
})

## Combine tables
tbl_combined <- tbl_stack(cox_tables) %>%
  bold_labels() %>%
  modify_caption("**Multivariate Cox Regression Analysis**")

## Export as image
ft <- as_flex_table(tbl_combined)
save_as_image(ft, path = "cox_results_temp.png", zoom = 2)

## Export PDF
pdf(output_pdf, width = 11, height = 8.5)
grid.newpage()
grid.raster(readPNG("cox_results_temp.png"))
dev.off()

message("Multivariate Cox regression completed successfully.")
message("Results saved to: ", output_pdf)
