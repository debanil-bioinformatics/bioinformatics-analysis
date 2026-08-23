library(tidyverse)
library(ggplot2)

# ── 1. LOAD ───────────────────────────────────────────────────────────────────

expr <- read_tsv("~/Documents/Cancers/brca_metabric/data_mrna_illumina_microarray.txt",
                 show_col_types = FALSE)

# Extract BUB1B and MAD2L1
bub1b <- expr %>% filter(Hugo_Symbol == "BUB1B") %>% select(-Hugo_Symbol, -Entrez_Gene_Id) %>% as.numeric()
mad2l1 <- expr %>% filter(Hugo_Symbol == "MAD2L1") %>% select(-Hugo_Symbol, -Entrez_Gene_Id) %>% as.numeric()

df <- tibble(
  BUB1B = bub1b,
  MAD2L1 = mad2l1
) %>%
  drop_na()

cat("Samples:", nrow(df), "\n")

# ── 2. CORRELATION ────────────────────────────────────────────────────────────

res <- cor.test(df$BUB1B, df$MAD2L1, method = "pearson")
r  <- res$estimate
p  <- res$p.value
r2 <- r^2

cat(sprintf("\nPearson r = %.3f\np-value   = %.2e\nR²        = %.3f\nn         = %d\n",
            r, p, r2, nrow(df)))

# ── 3. PLOT ───────────────────────────────────────────────────────────────────

label <- sprintf("r = %.3f\np = %.2e\nR² = %.3f\nn = %d", r, p, r2, nrow(df))

ggplot(df, aes(x = BUB1B, y = MAD2L1)) +
  geom_point(color = "#2c7bb6", size = 2, alpha = 0.5) +
  geom_smooth(method = "lm", color = "#d7191c", fill = "#d7191c",
              alpha = 0.15, linewidth = 0.9) +
  annotate("text", x = Inf, y = -Inf, label = label,
           hjust = 1.1, vjust = -0.5, size = 4, family = "mono") +
  labs(
    title = "BUB1B vs MAD2L1 — METABRIC",
    x = "BUB1B (microarray expression)",
    y = "MAD2L1 (microarray expression)",
    caption = "Source: METABRIC via cBioPortal"
  ) +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))

ggsave("metabric_bub1b_mad2l1.png", width = 7, height = 6, dpi = 300, bg = "white")
ggsave("metabric_bub1b_mad2l1.pdf", width = 7, height = 6)

cat("Done!\n")