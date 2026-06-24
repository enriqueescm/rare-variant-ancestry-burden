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
  filter(an >= 10000) %>%        # minimum coverage threshold
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

# ── 12. Ancestry-specific variants: lollipop plot ─────────────
# Define "ancestry-specific" as: variant observed in only 1 ancestry group
# with AC > 0, among all ancestry groups with sufficient coverage (AN > 1000)

# Step 1: pivot to wide format — one row per variant
ancestry_specific <- df_func %>%
  filter(an > 1000) %>%                         # only well-covered ancestries
  group_by(gene, variant_id, ancestry) %>%
  summarise(ac = sum(ac), .groups = "drop") %>%
  group_by(gene, variant_id) %>%
  summarise(
    n_ancestries_with_ac = sum(ac > 0),          # how many ancestries carry it
    dominant_ancestry    = ancestry[which.max(ac)],
    max_ac               = max(ac),
    .groups = "drop"
  ) %>%
  filter(n_ancestries_with_ac == 1)              # strictly ancestry-specific

cat(sprintf("\nAncestry-specific variants: %s\n", nrow(ancestry_specific)))

# Step 2: get mean AN per gene per ancestry for normalization
an_per_ancestry <- df_func %>%
  filter(an > 1000) %>%
  group_by(gene, ancestry) %>%
  summarise(mean_an = mean(an), .groups = "drop")

# Step 3: count specific variants per gene per ancestry, normalized
lollipop_data <- ancestry_specific %>%
  left_join(module_map, by = "gene") %>%
  count(gene, dominant_ancestry, module, name = "n_specific") %>%
  left_join(an_per_ancestry,
            by = c("gene", "dominant_ancestry" = "ancestry")) %>%
  filter(dominant_ancestry %in% ancestry_order) %>%
  mutate(
    # Normalize: specific variants per 10,000 chromosomes
    n_specific_norm = (n_specific / mean_an) * 1e4,
    ancestry_label  = ancestry_labels_full[dominant_ancestry],
    gene            = factor(gene, levels = gene_order),
    dominant_ancestry = factor(dominant_ancestry, levels = ancestry_order)
  )

cat("\nNormalized ancestry-specific variants (top 20):\n")
lollipop_data %>%
  arrange(desc(n_specific_norm)) %>%
  head(20) %>%
  select(gene, dominant_ancestry, n_specific, mean_an, n_specific_norm) %>%
  as.data.frame() %>%
  print()

# Step 4: plot normalized
p_lollipop <- lollipop_data %>%
  ggplot(aes(x = n_specific_norm, y = gene, color = dominant_ancestry)) +
  geom_segment(aes(x = 0, xend = n_specific_norm, y = gene, yend = gene),
               linewidth = 0.6, alpha = 0.6) +
  geom_point(size = 3.5, alpha = 0.9) +
  facet_grid(module ~ ., scales = "free_y", space = "free_y") +
  scale_color_manual(
    values = c(
      afr = "#e94560",
      ami = "#a8dadc",
      amr = "#f4a261",
      asj = "#9b5de5",
      eas = "#00b4d8",
      fin = "#06d6a0",
      mid = "#ffd166",
      nfe = "#457b9d",
      sas = "#f77f00"
    ),
    labels = ancestry_labels_full,
    name   = "Ancestry"
  ) +
  labs(
    title    = "Ancestry-Specific Rare Functional Variants (Normalized)",
    subtitle = "Variants present in exactly one ancestry | per 10,000 chromosomes | gnomAD v4",
    x        = "Ancestry-specific variants per 10,000 chromosomes",
    y        = NULL,
    caption  = "Normalized by mean AN per gene per ancestry to correct for sample size differences"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.y      = element_text(size = 9, face = "italic"),
    strip.text.y     = element_text(angle = 0, face = "bold", size = 8),
    strip.background = element_rect(fill = "#f0f0f0", color = NA),
    legend.position  = "right",
    plot.title       = element_text(face = "bold", size = 13),
    plot.subtitle    = element_text(size = 9, color = "grey40"),
    plot.caption     = element_text(size = 7, color = "grey60"),
    panel.grid.minor = element_blank()
  )

ggsave(
  "figures/02_ancestry_specific_lollipop.png",
  plot   = p_lollipop,
  width  = 10,
  height = 8,
  dpi    = 300
)

cat("\nLollipop plot saved: figures/02_ancestry_specific_lollipop.png\n")

# ── 13. Scatter: demographic history vs burden ────────────────
# Effective population size (Ne) estimates from published literature
# Sources: Tenesa et al. 2007, Browning et al. 2018, gnomAD flagship paper
# These are approximate historical Ne estimates (thousands)

ne_data <- tibble(
  ancestry = c("afr", "amr", "asj", "eas", "fin", "mid", "nfe", "sas"),
  ancestry_label = c("African", "Admixed American", "Ashkenazi Jewish",
                     "East Asian", "Finnish", "Middle Eastern",
                     "Non-Finnish Eur.", "South Asian"),
  # Historical Ne in thousands (approximate)
  # afr: large, pre-OOA; fin/asj: strong bottleneck
  ne_thousands = c(17.0, 5.0, 1.5, 7.0, 3.5, 4.0, 8.0, 9.0),
  # Approximate bottleneck strength (1=strong, 5=weak/none)
  bottleneck = c(5, 2, 1, 3, 1, 2, 3, 3)
)

# Calculate mean burden per ancestry across all genes
mean_burden <- burden %>%
  filter(ancestry %in% ne_data$ancestry) %>%
  group_by(ancestry) %>%
  summarise(
    mean_burden   = mean(burden_1e4, na.rm = TRUE),
    median_burden = median(burden_1e4, na.rm = TRUE),
    .groups = "drop"
  )

scatter_data <- ne_data %>%
  left_join(mean_burden, by = "ancestry")

cat("\nDemographic vs burden data:\n")
scatter_data %>%
  select(ancestry_label, ne_thousands, mean_burden) %>%
  arrange(desc(mean_burden)) %>%
  as.data.frame() %>%
  print()

# Plot
p_scatter <- scatter_data %>%
  ggplot(aes(x = ne_thousands, y = mean_burden)) +
  geom_smooth(method = "lm", se = TRUE,
              color = "#457b9d", fill = "#457b9d", alpha = 0.15,
              linewidth = 0.8) +
  geom_point(aes(fill = ancestry, size = bottleneck),
             shape = 21, color = "white", stroke = 0.5, alpha = 0.9) +
  geom_text(aes(label = ancestry_label),
            hjust = -0.15, vjust = 0.4, size = 3, color = "grey30") +
  scale_fill_manual(
    values = c(
      afr = "#e94560", amr = "#f4a261", asj = "#9b5de5",
      eas = "#00b4d8", fin = "#06d6a0", mid = "#ffd166",
      nfe = "#457b9d", sas = "#f77f00"
    ),
    guide = "none"
  ) +
  scale_size_continuous(
    range  = c(4, 10),
    name   = "Bottleneck\nstrength",
    breaks = c(1, 3, 5),
    labels = c("Strong", "Moderate", "Weak/None")
  ) +
  scale_x_continuous(
    name   = "Historical effective population size, Ne (thousands)",
    limits = c(0, 20)
  ) +
  labs(
    title    = "Demographic History Predicts Rare Variant Burden",
    subtitle = "Mean rare functional variant burden vs. historical Ne across ancestries",
    y        = "Mean burden (per 10K chromosomes)",
    caption  = "Ne estimates from Tenesa et al. 2007 and gnomAD flagship paper.\nPoint size reflects bottleneck strength."
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title    = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 9, color = "grey40"),
    plot.caption  = element_text(size = 7, color = "grey60"),
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

ggsave(
  "figures/03_demographic_burden_scatter.png",
  plot   = p_scatter,
  width  = 10,
  height = 7,
  dpi    = 300
)

cat("\nScatter plot saved: figures/03_demographic_scatter.png\n")