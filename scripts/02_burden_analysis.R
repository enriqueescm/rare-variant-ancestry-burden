# 02_burden_analysis.R
# Ancestry-stratified rare variant burden analysis
# gnomAD v4 | 15 candidate genes | 3 biological modules

# ── 1. Libraries ──────────────────────────────────────────────
library(tidyverse)
library(data.table)

# ── 2. Load data ──────────────────────────────────────────────
cat("Loading data...\n")

df <- fread(
  "data/raw/all_genes_gnomad_v4.tsv",
  colClasses = list(character = c("lof", "lof_filter", "consequence"))
)

cat(sprintf("Total rows loaded: %s\n", format(nrow(df), big.mark = ",")))
cat(sprintf("Genes: %s\n", paste(unique(df$gene), collapse = ", ")))
cat(sprintf("Ancestries: %s\n", paste(unique(df$ancestry), collapse = ", ")))

# ── 3. Filter: exomes only, per-ancestry rows, exclude ALL ────
cat("\nFiltering...\n")

df_clean <- df %>%
  filter(
    data_type == "exome",
    ancestry  != "ALL",
    ancestry  != "remaining",
    an        > 0
  )

# ── 4. Filter: rare variants only (AF < 1%) ───────────────────
df_rare <- df_clean %>%
  filter(af < 0.01)

cat(sprintf("Rare variants (AF < 1%%): %s rows\n", format(nrow(df_rare), big.mark = ",")))

# ── 5. Filter: functional variants only ──────────────────────
df_func <- df_rare %>%
  filter(
    lof == "HC" |
      consequence %in% c(
        "missense_variant",
        "stop_gained",
        "frameshift_variant",
        "splice_donor_variant",
        "splice_acceptor_variant",
        "start_lost",
        "stop_lost"
      )
  )

cat(sprintf("Functional rare variants: %s rows\n", format(nrow(df_func), big.mark = ",")))
cat(sprintf("Consequence types: %s\n", paste(unique(df_func$consequence), collapse = ", ")))

# ── 6. Sanity check: variants per gene ───────────────────────
cat("\nVariants per gene (functional, rare):\n")

df_func %>%
  distinct(gene, variant_id) %>%
  count(gene, name = "n_variants") %>%
  arrange(desc(n_variants)) %>%
  as.data.frame() %>%
  print()

  # ── 7. Calculate burden per gene per ancestry ─────────────────
# Burden = sum of rare functional allele counts (AC)
#          normalized by total allele number (AN)
# This gives us: expected rare functional variants per chromosome
# Units: variants per 10,000 chromosomes (x10^4 for readability)

burden <- df_func %>%
  group_by(gene, ancestry) %>%
  summarise(
    total_ac    = sum(ac, na.rm = TRUE),    # total rare functional alleles
    mean_an     = mean(an, na.rm = TRUE),   # mean coverage (proxy for sample size)
    n_variants  = n_distinct(variant_id),   # unique variants contributing
    burden_raw  = total_ac / mean_an,       # raw burden
    burden_1e4  = burden_raw * 1e4,         # scaled per 10,000 chromosomes
    .groups = "drop"
  )

cat("\nBurden summary (top 20 gene x ancestry combinations):\n")
burden %>%
  arrange(desc(burden_1e4)) %>%
  head(20) %>%
  as.data.frame() %>%
  print()

# ── 8. Add module annotation ──────────────────────────────────
module_map <- tibble(
  gene = c("IKBKB", "IFIH1", "TLR1", "TLR6", "TLR10", "ACKR1", "CD36",
           "STAT4", "IRF5", "TNFAIP3",
           "RUNX2", "COL1A1", "SP7", "LRP5", "VDR"),
  module = c(rep("A: Innate immunity", 7),
             rep("B: Autoimmunity", 3),
             rep("C: Skeletal", 5))
)

burden <- burden %>%
  left_join(module_map, by = "gene")

# ── 9. Ancestry labels (readable names) ──────────────────────
ancestry_labels <- c(
  afr = "African",
  ami = "Amish",
  amr = "Admixed American",
  asj = "Ashkenazi Jewish",
  eas = "East Asian",
  fin = "Finnish",
  mid = "Middle Eastern",
  nfe = "Non-Finnish European",
  sas = "South Asian"
)

burden <- burden %>%
  mutate(ancestry_label = ancestry_labels[ancestry])

# ── 10. Save processed data ───────────────────────────────────
write_tsv(burden, "data/processed/burden_by_gene_ancestry.tsv")
cat(sprintf("\nBurden table saved: %s gene x ancestry combinations\n",
            nrow(burden)))

# ── 11. Heatmap: burden × gene × ancestry ────────────────────
library(ggplot2)

# Check ancestries present
cat("\nAncestries in burden table:\n")
print(unique(burden$ancestry))

# Order genes by burden in Africans (descending) within each module
gene_order <- burden %>%
  filter(ancestry == "afr") %>%
  arrange(module, desc(burden_1e4)) %>%
  pull(gene)

ancestry_order <- c("afr", "ami", "amr", "asj", "eas", "fin", "mid", "nfe", "sas")

ancestry_labels_full <- c(
  afr = "African",
  ami = "Amish",
  amr = "Admixed American",
  asj = "Ashkenazi Jewish",
  eas = "East Asian",
  fin = "Finnish",
  mid = "Middle Eastern",
  nfe = "Non-Finnish Eur.",
  sas = "South Asian"
)

p_heatmap <- burden %>%
  filter(ancestry %in% ancestry_order) %>%
  mutate(
    gene          = factor(gene, levels = gene_order),
    ancestry      = factor(ancestry, levels = ancestry_order),
    burden_capped = pmin(burden_1e4, quantile(burden_1e4, 0.95)),
    label_text    = ifelse(gene == "CD36" & ancestry == "mid",
                           paste0(round(burden_1e4, 0), "*"),
                           as.character(round(burden_1e4, 0)))
  ) %>%
  ggplot(aes(x = ancestry, y = gene, fill = burden_capped)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = label_text),
            size = 2.5, color = "white", fontface = "bold") +
  facet_grid(module ~ ., scales = "free_y", space = "free_y") +
  scale_fill_gradientn(
    colors = c("#1a1a2e", "#16213e", "#0f3460", "#533483", "#e94560"),
    name   = "Burden\n(per 10K chr)"
  ) +
  scale_x_discrete(labels = ancestry_labels_full) +
  labs(
    title    = "Rare Functional Variant Burden Across Human Populations",
    subtitle = "gnomAD v4 | MAF < 1% | LoF HC + missense | 15 candidate genes",
    x        = NULL,
    y        = NULL,
    caption  = "Burden = Σ(AC) / mean(AN) × 10,000  |  * small sample size (AN < 10,000)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x      = element_text(angle = 35, hjust = 1, size = 9),
    axis.text.y      = element_text(size = 9, face = "italic"),
    strip.text.y     = element_text(angle = 0, face = "bold", size = 8),
    strip.background = element_rect(fill = "#f0f0f0", color = NA),
    legend.position  = "right",
    plot.title       = element_text(face = "bold", size = 13),
    plot.subtitle    = element_text(size = 9, color = "grey40"),
    plot.caption     = element_text(size = 7, color = "grey60"),
    panel.grid       = element_blank()
  )

ggsave(
  "figures/01_burden_heatmap.png",
  plot   = p_heatmap,
  width  = 11,
  height = 8,
  dpi    = 300
)

cat("\nHeatmap saved: figures/01_burden_heatmap.png\n")