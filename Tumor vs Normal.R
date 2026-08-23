## --------------------------------------------------
## CUEDC2 expression analysis
## Tumour vs Normal (overall) and HR stratified
## + bracketed p-value + statistics export
## --------------------------------------------------

## 1. Set working directory
setwd("~/Documents/TCGA Analysis")

## 2. Load libraries
library(dplyr)
library(ggplot2)
library(ggpubr)

## 3. Read annotated dataset
tcga <- read.csv("TCGA_BRCA_CUEDC2_annotated.csv", stringsAsFactors = FALSE)

## --------------------------------------------------
## A. Tumour vs Normal – whole cohort
## --------------------------------------------------

tcga_tn <- tcga %>%
  filter(
    sample_type %in% c("Tumour", "Normal"),
    !is.na(CUEDC2)
  )

## Statistical test
wilcox_tn <- wilcox.test(CUEDC2 ~ sample_type, data = tcga_tn)

## Medians
median_tumour <- median(tcga_tn$CUEDC2[tcga_tn$sample_type == "Tumour"])
median_normal <- median(tcga_tn$CUEDC2[tcga_tn$sample_type == "Normal"])
median_diff_tn <- median_tumour - median_normal

## Plot (with bracket)
p_tn <- ggplot(tcga_tn, aes(x = sample_type, y = CUEDC2, fill = sample_type)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  geom_jitter(width = 0.15, size = 0.35, alpha = 0.4) +
  stat_compare_means(
    comparisons = list(c("Tumour", "Normal")),
    method = "wilcox.test",
    label = "p.signif",
    tip.length = 0.02
  ) +
  theme_classic(base_size = 14) +
  labs(
    x = "",
    y = "CUEDC2 expression (log2(norm_count + 1))",
    title = "CUEDC2 expression: Tumour vs Normal"
  ) +
  theme(legend.position = "none")

print(p_tn)

ggsave("CUEDC2_Tumour_vs_Normal.pdf", p_tn, width = 5, height = 5)
ggsave("CUEDC2_Tumour_vs_Normal.tiff", p_tn,
       width = 5, height = 5, units = "in",
       dpi = 600, compression = "lzw")

## --------------------------------------------------
## B. CUEDC2 expression by HR status (tumour only)
## --------------------------------------------------

tcga_hr <- tcga %>%
  filter(
    sample_type == "Tumour",
    HR_status %in% c("HR+", "HR-"),
    !is.na(CUEDC2)
  )

## Statistical test
wilcox_hr <- wilcox.test(CUEDC2 ~ HR_status, data = tcga_hr)

## Medians
median_hr_pos <- median(tcga_hr$CUEDC2[tcga_hr$HR_status == "HR+"])
median_hr_neg <- median(tcga_hr$CUEDC2[tcga_hr$HR_status == "HR-"])
median_diff_hr <- median_hr_neg - median_hr_pos

## Plot (with bracket)
p_hr <- ggplot(tcga_hr, aes(x = HR_status, y = CUEDC2, fill = HR_status)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  geom_jitter(width = 0.15, size = 0.35, alpha = 0.4) +
  stat_compare_means(
    comparisons = list(c("HR+", "HR-")),
    method = "wilcox.test",
    label = "p.signif",
    tip.length = 0.02
  ) +
  theme_classic(base_size = 14) +
  labs(
    x = "",
    y = "CUEDC2 expression (log2(norm_count + 1))",
    title = "CUEDC2 expression in tumours: HR+ vs HR−"
  ) +
  theme(legend.position = "none")

print(p_hr)

ggsave("CUEDC2_HRpos_vs_HRneg_Tumours.pdf", p_hr, width = 5, height = 5)
ggsave("CUEDC2_HRpos_vs_HRneg_Tumours.tiff", p_hr,
       width = 5, height = 5, units = "in",
       dpi = 600, compression = "lzw")

## --------------------------------------------------
## C. Statistical summary table (combined)
## --------------------------------------------------

stats_summary <- data.frame(
  Comparison = c("Tumour vs Normal", "HR- vs HR+ (Tumours)"),
  Test = c("Wilcoxon rank-sum test", "Wilcoxon rank-sum test"),
  Median_Group1 = c(median_tumour, median_hr_neg),
  Median_Group2 = c(median_normal, median_hr_pos),
  Median_Difference = c(median_diff_tn, median_diff_hr),
  P_value = c(wilcox_tn$p.value, wilcox_hr$p.value)
)

write.csv(
  stats_summary,
  file = "CUEDC2_expression_statistics_summary.csv",
  row.names = FALSE
)

## --------------------------------------------------
## End of analysis
## --------------------------------------------------
