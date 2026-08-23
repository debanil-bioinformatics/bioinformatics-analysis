library(tidyverse)
library(ggplot2)

# ── 1. LOAD FILES ─────────────────────────────────────────────────────────────

expr <- read_csv("Expression_(Short-read)_Public_26Q1_subsetted.csv",
                 show_col_types = FALSE) %>%
  rename(ModelID = 1)

meta <- read_csv("Model.csv", show_col_types = FALSE)

# ── 2. FILTER BREAST CANCER LINES ────────────────────────────────────────────

bc_ids <- meta %>%
  filter(OncotreeLineage == "Breast") %>%
  pull(ModelID)

cat("Breast cancer lines in metadata:", length(bc_ids), "\n")

# ── 3. EXTRACT BUB1B AND MAD2L1 ──────────────────────────────────────────────

df <- expr %>%
  filter(ModelID %in% bc_ids) %>%
  select(ModelID, BUB1B, MAD2L1) %>%
  drop_na() %>%
  left_join(meta %>% select(ModelID, CellLineName, OncotreeSubtype), by = "ModelID")

cat("Cell lines used for correlation:", nrow(df), "\n")

# ── 4. CELL LINE LIST ─────────────────────────────────────────────────────────

cat("\nCell lines included:\n")
print(df %>% select(CellLineName, OncotreeSubtype, BUB1B, MAD2L1), n = Inf)

# ── 5. CORRELATION ────────────────────────────────────────────────────────────

res <- cor.test(df$BUB1B, df$MAD2L1, method = "pearson")
r   <- res$estimate
p   <- res$p.value
r2  <- r^2

cat(sprintf("\nPearson r = %.3f\np-value   = %.2e\nR²        = %.3f\nn         = %d\n",
            r, p, r2, nrow(df)))

# ── 6. PLOT ───────────────────────────────────────────────────────────────────

label <- sprintf("r = %.3f\np = %.2e\nR² = %.3f\nn = %d", r, p, r2, nrow(df))

ggplot(df, aes(x = BUB1B, y = MAD2L1)) +
  geom_point(color = "#2c7bb6", size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", color = "#d7191c", fill = "#d7191c",
              alpha = 0.15, linewidth = 0.9) +
  annotate("text", x = Inf, y = -Inf, label = label,
           hjust = 1.1, vjust = -0.5, size = 4, family = "mono") +
  labs(
    title = "BUB1B vs MAD2L1 — CCLE Breast Cancer Cell Lines",
    x = "BUB1B (log₂ TPM+1)",
    y = "MAD2L1 (log₂ TPM+1)",
    caption = "Source: DepMap 26Q1"
  ) +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))

ggsave("ccle_bub1b_mad2l1.png", width = 7, height = 6, dpi = 300, bg = "white")
ggsave("ccle_bub1b_mad2l1.pdf", width = 7, height = 6)

# ── 7. SAVE CELL LINE LIST ────────────────────────────────────────────────────

write_csv(df %>% select(CellLineName, OncotreeSubtype, BUB1B, MAD2L1),
          "ccle_celllines_used.csv")

cat("\nDone! Plot and cell line list saved.\n")