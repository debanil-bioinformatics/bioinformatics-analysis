library(tidyverse)
library(ggplot2)

# ── 1. LOAD ───────────────────────────────────────────────────────────────────

df <- read_csv("TCGA BRCA MAD2L1 BUB1B ESR1.csv", show_col_types = FALSE) %>%
  drop_na(ESR1, BUB1B)

cat("Samples:", nrow(df), "\n")
cat("Sample types:\n")
print(table(df$sample_type))

# ── 2. FILTER PRIMARY TUMOR ONLY ──────────────────────────────────────────────

df <- df %>% filter(sample_type == "Primary Tumor")
cat("Primary tumor samples:", nrow(df), "\n")

# ── 3. CORRELATION ────────────────────────────────────────────────────────────

res <- cor.test(df$ESR1, df$BUB1B, method = "pearson")
r  <- res$estimate
p  <- res$p.value
r2 <- r^2

cat(sprintf("\nPearson r = %.3f\np-value   = %.2e\nR²        = %.3f\nn         = %d\n",
            r, p, r2, nrow(df)))

# ── 4. PLOT ───────────────────────────────────────────────────────────────────

label <- sprintf("r = %.3f\np = %.2e\nR² = %.3f\nn = %d", r, p, r2, nrow(df))

ggplot(df, aes(x = ESR1, y = BUB1B)) +
  geom_point(color = "#2c7bb6", size = 2, alpha = 0.5) +
  geom_smooth(method = "lm", color = "#d7191c", fill = "#d7191c",
              alpha = 0.15, linewidth = 0.9) +
  annotate("text", x = Inf, y = -Inf, label = label,
           hjust = 1.1, vjust = -0.5, size = 4, family = "mono") +
  labs(
    title = "ESR1 vs BUB1B — TCGA BRCA",
    x = "ESR1 (log₂ normalized expression)",
    y = "BUB1B (log₂ normalized expression)",
    caption = "Source: TCGA BRCA via UCSC Xena"
  ) +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))

ggsave("tcga_ESR1_BUB1B.png", width = 7, height = 6, dpi = 300, bg = "white")
ggsave("tcga_ESR1_BUB1B.pdf", width = 7, height = 6)

cat("Done!\n")