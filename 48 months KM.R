## --------------------------------------------------
## HR− only survival analysis
## CUEDC2 stratified within HR− cohort (median-based)
## KM OS & PFI (MONTHS)
## 48-MONTH ADMINISTRATIVE CENSORING (EVENT-PRESERVING)
## KM plots WITHOUT risk table
## Risk table INCLUDED in KM statistics CSV
## --------------------------------------------------

## 1. Set working directory
setwd("~/Documents/TCGA Analysis")

## 2. Load libraries
library(dplyr)
library(survival)
library(survminer)
library(tidyr)

## --------------------------------------------------
## A. Read annotated dataset
## --------------------------------------------------

tcga <- read.csv("TCGA HR NEG.csv", stringsAsFactors = FALSE)

tcga_hrneg <- tcga %>%
  filter(
    sample_type == "Tumour",
    HR_status == "HR-",
    !is.na(CUEDC2)
  )

## --------------------------------------------------
## B. Define CUEDC2 High / Low (HR− median)
## --------------------------------------------------

median_cuedc2_hrneg <- median(tcga_hrneg$CUEDC2, na.rm = TRUE)

tcga_hrneg <- tcga_hrneg %>%
  mutate(
    CUEDC2_group = ifelse(CUEDC2 >= median_cuedc2_hrneg, "High", "Low"),
    CUEDC2_group = factor(CUEDC2_group, levels = c("High", "Low"))
  )

## --------------------------------------------------
## Global administrative censoring at 48 months
## --------------------------------------------------

MAX_FU <- 48  # months

time_points <- c(0, 12, 24, 36, 48)

## --------------------------------------------------
## C. Overall Survival (OS) — 48 month follow-up capped
## --------------------------------------------------

tcga_os <- tcga_hrneg %>%
  filter(!is.na(OS), !is.na(OS.time)) %>%
  mutate(
    OS_months_raw = OS.time / 30.44,
    ## TRUE administrative censoring
    OS_event_48   = ifelse(OS_months_raw > MAX_FU, 0, OS),
    OS_months_48  = pmin(OS_months_raw, MAX_FU)
  )

fit_os <- survfit(Surv(OS_months_48, OS_event_48) ~ CUEDC2_group, data = tcga_os)

## KM plot (NO risk table)
p_km_os <- ggsurvplot(
  fit_os,
  data = tcga_os,
  risk.table = FALSE,
  pval = TRUE,
  conf.int = FALSE,
  xlab = "Time (months)",
  ylab = "Overall survival probability",
  legend.title = "CUEDC2",
  legend.labs = c("High", "Low"),
  palette = c("red", "black"),
  ggtheme = theme_classic(base_size = 14),
  xlim = c(0, 48)
)

print(p_km_os)

ggsave("KM48_OS_HRneg_CUEDC2_HRmedian.pdf",
       p_km_os$plot, width = 6, height = 5)

ggsave("KM48_OS_HRneg_CUEDC2_HRmedian.tiff",
       p_km_os$plot, width = 6, height = 5,
       units = "in", dpi = 600, compression = "lzw")

## Log-rank test (OS, 48m capped)
lr_os <- survdiff(Surv(OS_months_48, OS_event_48) ~ CUEDC2_group, data = tcga_os)
p_os <- 1 - pchisq(lr_os$chisq, df = 1)

## Risk table extraction (OS)
os_risk <- summary(fit_os, times = time_points)

risk_os_df <- data.frame(
  Endpoint = "OS",
  Group = gsub("CUEDC2_group=", "", os_risk$strata),
  Time_months = os_risk$time,
  N_at_risk = os_risk$n.risk
) %>%
  pivot_wider(
    names_from = Time_months,
    values_from = N_at_risk,
    names_prefix = "Month_"
  )

## --------------------------------------------------
## D. Progression-Free Interval (PFI) — 48 month follow-up capped
## --------------------------------------------------

tcga_pfi <- tcga_hrneg %>%
  filter(!is.na(PFI), !is.na(PFI.time)) %>%
  mutate(
    PFI_months_raw = PFI.time / 30.44,
    ## TRUE administrative censoring
    PFI_event_48   = ifelse(PFI_months_raw > MAX_FU, 0, PFI),
    PFI_months_48  = pmin(PFI_months_raw, MAX_FU)
  )

fit_pfi <- survfit(Surv(PFI_months_48, PFI_event_48) ~ CUEDC2_group, data = tcga_pfi)

## KM plot (NO risk table)
p_km_pfi <- ggsurvplot(
  fit_pfi,
  data = tcga_pfi,
  risk.table = FALSE,
  pval = TRUE,
  conf.int = FALSE,
  xlab = "Time (months)",
  ylab = "Progression-free interval probability",
  legend.title = "CUEDC2",
  legend.labs = c("High", "Low"),
  palette = c("red", "black"),
  ggtheme = theme_classic(base_size = 14),
  xlim = c(0, 48)
)

print(p_km_pfi)

ggsave("KM48_PFI_HRneg_CUEDC2_HRmedian.pdf",
       p_km_pfi$plot, width = 6, height = 5)

ggsave("KM48_PFI_HRneg_CUEDC2_HRmedian.tiff",
       p_km_pfi$plot, width = 6, height = 5,
       units = "in", dpi = 600, compression = "lzw")

## Log-rank test (PFI, 48m capped)
lr_pfi <- survdiff(Surv(PFI_months_48, PFI_event_48) ~ CUEDC2_group, data = tcga_pfi)
p_pfi <- 1 - pchisq(lr_pfi$chisq, df = 1)

## Risk table extraction (PFI)
pfi_risk <- summary(fit_pfi, times = time_points)

risk_pfi_df <- data.frame(
  Endpoint = "PFI",
  Group = gsub("CUEDC2_group=", "", pfi_risk$strata),
  Time_months = pfi_risk$time,
  N_at_risk = pfi_risk$n.risk
) %>%
  pivot_wider(
    names_from = Time_months,
    values_from = N_at_risk,
    names_prefix = "Month_"
  )

## --------------------------------------------------
## E. Combine KM statistics + risk tables (SINGLE CSV)
## --------------------------------------------------

km_stats <- bind_rows(risk_os_df, risk_pfi_df) %>%
  mutate(
    Test = "Log-rank test",
    Cutoff = "Median within HR− cohort",
    Followup = "Administratively censored at 48 months",
    P_value = ifelse(Endpoint == "OS", p_os, p_pfi)
  )

write.csv(
  km_stats,
  file = "HRneg_CUEDC2_KM48_statistics_withRiskTable.csv",
  row.names = FALSE
)

## --------------------------------------------------
## End of script
## --------------------------------------------------
