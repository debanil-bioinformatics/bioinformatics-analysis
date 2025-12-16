############################################################
# Author: Debanil Dhar
# Year: 2025
#
# Description:
# Generalized Kaplan–Meier survival analysis and Cox
# proportional hazards modeling for any cancer cohort.
# Supports user-defined biomarkers, survival endpoints,
# optional time censoring, and publication-ready output.
#
# Input:
# - CSV file with survival time, event status, and a
#   grouping (biomarker) variable
#
# Output:
# - Kaplan–Meier survival plot (PDF)
# - Supplementary PDF with Cox statistics and risk table
############################################################

## Load required libraries
suppressPackageStartupMessages({
  library(survival)
  library(survminer)
  library(dplyr)
  library(ggplot2)
})

## ---------------- USER PARAMETERS ----------------

# Input
input_file <- "survival_data.csv"

# Biomarker / grouping variable
group_col    <- "Biomarker_Status"   # e.g., High / Low
group_levels <- c("Low", "High")      # reference first

# Survival endpoint
time_col  <- "OS_time"    # survival time column
event_col <- "OS_event"   # 1 = event, 0 = censored

# Optional censoring (set to NA to disable)
max_followup <- 200       # e.g., months

# Plot labels
plot_title <- "Overall Survival by Biomarker Status"
x_label    <- "Time (Months)"
y_label    <- "Survival Probability"

# Output files
km_pdf      <- "kaplan_meier_plot.pdf"
stats_pdf   <- "survival_statistics.pdf"

## -------------------------------------------------

## Load data
data <- read.csv(input_file, stringsAsFactors = FALSE)

## Sanity checks
stopifnot(all(c(group_col, time_col, event_col) %in% colnames(data)))

## Factorize grouping variable
data[[group_col]] <- factor(data[[group_col]], levels = group_levels)

## Apply optional censoring
if (!is.na(max_followup)) {
  data$survival_time <- pmin(data[[time_col]], max_followup)
} else {
  data$survival_time <- data[[time_col]]
}

data$event_status <- data[[event_col]]

## Create survival object
surv_obj <- Surv(time = data$survival_time, event = data$event_status)

## Kaplan–Meier fit
surv_fit <- survfit(surv_obj ~ data[[group_col]], data = data)

## Cox proportional hazards model
cox_model <- coxph(surv_obj ~ data[[group_col]], data = data)
cox_summary <- summary(cox_model)

## Extract HR and p-value
coef_name <- paste0("data[[group_col]]", group_levels[2])
hr <- round(cox_summary$coefficients[1, "exp(coef)"], 2)
p_value <- format.pval(cox_summary$coefficients[1, "Pr(>|z|)"], digits = 3)

## Create Kaplan–Meier plot
km_plot <- ggsurvplot(
  surv_fit,
  data = data,
  conf.int = FALSE,
  risk.table = FALSE,
  pval = FALSE,
  palette = c("black", "red"),
  title = plot_title,
  xlab = x_label,
  ylab = y_label,
  legend.title = group_col,
  legend.labs = group_levels,
  xlim = if (!is.na(max_followup)) c(0, max_followup) else NULL
)

## Annotate HR and p-value
km_plot$plot <- km_plot$plot +
  annotate(
    "text",
    x = max(data$survival_time, na.rm = TRUE) * 0.7,
    y = 0.95,
    label = paste("HR (", group_levels[2], " vs ", group_levels[1], "): ", hr,
                  "\np-value: ", p_value, sep = ""),
    hjust = 0,
    size = 4
  )

## Save KM plot
ggsave(km_pdf, km_plot$plot, width = 8, height = 6)

## Risk table at fixed times
times_of_interest <- c(0, max_followup / 2, max_followup)
surv_summary <- summary(surv_fit, times = times_of_interest)

risk_table <- data.frame(
  Time = times_of_interest,
  Low  = surv_summary$n.risk[surv_summary$strata == paste0("data[[group_col]]=", group_levels[1])],
  High = surv_summary$n.risk[surv_summary$strata == paste0("data[[group_col]]=", group_levels[2])]
)

## Cox model text output
cox_text <- paste(capture.output(print(cox_summary)), collapse = "\n")

## Save supplementary PDF
pdf(stats_pdf, width = 8, height = 8)
par(mar = c(2, 2, 2, 2))
plot.new()

text(0.05, 0.95,
     paste("Risk Table:", paste(times_of_interest, "months", collapse = ", ")),
     adj = c(0, 1), family = "mono", cex = 1.1)

text(0.05, 0.90,
     paste(group_levels[2], ":", paste(risk_table$High, collapse = "  ")),
     adj = c(0, 1), family = "mono", cex = 1.1, col = "red")

text(0.05, 0.85,
     paste(group_levels[1], " :", paste(risk_table$Low, collapse = "  ")),
     adj = c(0, 1), family = "mono", cex = 1.1, col = "black")

text(0.05, 0.75, cox_text,
     adj = c(0, 1), family = "mono", cex = 0.75)

dev.off()

message("Kaplan–Meier and Cox analysis completed successfully.")
