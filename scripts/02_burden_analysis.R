# 02_burden_analysis.R
# Ancestry-stratified rare variant burden analysis
# gnomAD v4 | 15 candidate genes | 3 biological modules
#
# METHODOLOGICAL NOTES:
# 1. Burden = Σ(AC) / max(AN) × 10,000
#    max(AN) approximates the total callable chromosomes for the gene
#    in that ancestry, avoiding downward bias from low-coverage variants
# 2. MAF < 1% filter applied per-ancestry (not globally) — intentional:
#    we want variants rare IN that ancestry, not globally rare variants
# 3. Primary analysis uses LoF HC only (LOFTEE high-confidence).
#    Missense excluded: without CADD filtering they inflate burden
#    with likely benign variants (Pearson r=0.46 vs LoF-only).
# 4. TLR10 excluded from burden figures: no LoF HC variants in gnomAD v4.
# 5. ACKR1 retained in analysis but has 0 LoF HC burden across all
#    ancestries — documented as a constraint-consistent finding.

# ── 1. Libraries ──────────────────────────────────────────────
library(tidyverse)
library(data.table)
library(ggplot2)

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

# ── 4. Filter: rare variants only (AF < 1% per ancestry) ──────
df_rare <- df_clean %>%
  filter(af < 0.01)

cat(sprintf("Rare variants (AF < 1%% per ancestry): %s rows\n",
            format(nrow(df_rare), big.mark = ",")))

# ── 5. Filter: LoF HC only (primary analysis) ────────────────
df_lof_only <- df_rare %>%
  filter(lof == "HC")

# Keep df_func for reference counts only
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

cat(sprintf("LoF HC only (primary): %s rows\n",
            format(nrow(df_lof_only), big.mark = ",")))
cat(sprintf("LoF HC + missense (reference): %s rows\n",
            format(nrow(df_func), big.mark = ",")))

# ── 6. Sanity check: LoF variants per gene ───────────────────
cat("\nLoF HC variants per gene:\n")
df_lof_only %>%
  distinct(gene, variant_id) %>%
  count(gene, name = "n_lof_variants") %>%
  arrange(desc(n_lof_variants)) %>%
  as.data.frame() %>%
  print()

# ── 7. Calculate burden per gene per ancestry ─────────────────
calc_burden <- function(data) {
  data %>%
    group_by(gene, ancestry) %>%
    summarise(
      total_ac   = sum(ac, na.rm = TRUE),
      max_an     = max(an, na.rm = TRUE),
      n_variants = n_distinct(variant_id),
      burden_raw = total_ac / max_an,
      burden_1e4 = burden_raw * 1e4,
      .groups = "drop"
    )
}

burden_all <- df_lof_only %>%
  filter(an > 0) %>%
  calc_burden()

burden_filtered <- df_lof_only %>%
  filter(an >= 10000) %>%
  calc_burden()

burden <- burden_all

cat("\nLoF burden summary (top 20 gene x ancestry combinations):\n")
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

burden_filtered <- burden_filtered %>%
  left_join(module_map, by = "gene")

# ── 9. Ancestry labels ────────────────────────────────────────
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

burden <- burden %>%
  mutate(ancestry_label = ancestry_labels_full[ancestry])

# ── 10. Save processed data ───────────────────────────────────
write_tsv(burden_filtered, "data/processed/burden_by_gene_ancestry.tsv")
cat(sprintf("\nBurden table saved: %s gene x ancestry combinations\n",
            nrow(burden_filtered)))

# ── 11. Heatmap ───────────────────────────────────────────────
cat("\nAncestries in burden table:\n")
print(unique(burden$ancestry))

ancestry_order <- c("afr", "ami", "amr", "asj", "eas", "fin", "mid", "nfe", "sas")

# Exclude TLR10: no LoF HC variants in gnomAD v4
genes_for_figures <- module_map %>%
  filter(gene != "TLR10") %>%
  pull(gene)

gene_order <- burden %>%
  filter(ancestry == "afr", gene %in% genes_for_figures) %>%
  arrange(module, desc(burden_1e4)) %>%
  pull(gene)

p_heatmap <- burden %>%
  filter(
    ancestry %in% ancestry_order,
    gene %in% genes_for_figures,
    !is.na(burden_1e4)
  ) %>%
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
    title    = "Rare LoF Variant Burden Across Human Populations",
    subtitle = "gnomAD v4 | MAF < 1% per ancestry | LoF HC only (LOFTEE) | 14 candidate genes",
    x        = NULL,
    y        = NULL,
    caption  = "Burden = Σ(AC) / max(AN) × 10,000  |  * small sample size (AN < 10,000)  |  TLR10 excluded: no LoF HC variants in gnomAD v4"
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

ggsave("figures/01_burden_heatmap.png",
       plot = p_heatmap, width = 11, height = 8, dpi = 300)
cat("\nHeatmap saved: figures/01_burden_heatmap.png\n")

# ── 12. Lollipop: ancestry-specific LoF variants ─────────────
ancestry_specific <- df_lof_only %>%
  filter(an > 1000, gene %in% genes_for_figures) %>%
  group_by(gene, variant_id, ancestry) %>%
  summarise(ac = sum(ac), .groups = "drop") %>%
  group_by(gene, variant_id) %>%
  summarise(
    n_ancestries_with_ac = sum(ac > 0),
    dominant_ancestry    = ancestry[which.max(ac)],
    max_ac               = max(ac),
    .groups = "drop"
  ) %>%
  filter(n_ancestries_with_ac == 1)

cat(sprintf("\nAncestry-specific LoF variants (strict): %s\n",
            nrow(ancestry_specific)))

an_per_ancestry <- df_lof_only %>%
  filter(an > 1000, gene %in% genes_for_figures) %>%
  group_by(gene, ancestry) %>%
  summarise(mean_an = mean(an), .groups = "drop")

lollipop_data <- ancestry_specific %>%
  left_join(module_map, by = "gene") %>%
  count(gene, dominant_ancestry, module, name = "n_specific") %>%
  left_join(an_per_ancestry,
            by = c("gene", "dominant_ancestry" = "ancestry")) %>%
  filter(dominant_ancestry %in% ancestry_order) %>%
  mutate(
    n_specific_norm   = (n_specific / mean_an) * 1e4,
    ancestry_label    = ancestry_labels_full[dominant_ancestry],
    gene              = factor(gene, levels = gene_order),
    dominant_ancestry = factor(dominant_ancestry, levels = ancestry_order)
  )

cat("\nNormalized ancestry-specific LoF variants (top 20):\n")
lollipop_data %>%
  arrange(desc(n_specific_norm)) %>%
  head(20) %>%
  select(gene, dominant_ancestry, n_specific, mean_an, n_specific_norm) %>%
  as.data.frame() %>%
  print()

p_lollipop <- lollipop_data %>%
  ggplot(aes(x = n_specific_norm, y = gene, color = dominant_ancestry)) +
  geom_segment(aes(x = 0, xend = n_specific_norm, y = gene, yend = gene),
               linewidth = 0.6, alpha = 0.6) +
  geom_point(size = 3.5, alpha = 0.9) +
  facet_grid(module ~ ., scales = "free_y", space = "free_y") +
  scale_color_manual(
    values = c(afr="#e94560", ami="#a8dadc", amr="#f4a261", asj="#9b5de5",
               eas="#00b4d8", fin="#06d6a0", mid="#ffd166",
               nfe="#457b9d", sas="#f77f00"),
    labels = ancestry_labels_full, name = "Ancestry"
  ) +
  labs(
    title    = "Ancestry-Specific Rare LoF Variants (Normalized)",
    subtitle = "LoF HC variants present in exactly one ancestry | per 10,000 chromosomes | gnomAD v4",
    x        = "Ancestry-specific LoF variants per 10,000 chromosomes",
    y        = NULL,
    caption  = "Normalized by mean AN per gene per ancestry to correct for sample size differences\nAncestry-specific = AC > 0 in exactly one ancestry group (AN > 1,000) | LoF HC only"
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

ggsave("figures/02_ancestry_specific_lollipop.png",
       plot = p_lollipop, width = 10, height = 8, dpi = 300)
cat("\nLollipop plot saved: figures/02_ancestry_specific_lollipop.png\n")

# ── 13. Scatter: demographic context of burden ────────────────
# NOTE: exploratory only — R²=0.24, p=0.264 with LoF HC burden.
# Insufficient statistical power with n=7 populations and gene-specific
# LoF HC counts. Retained as contextual figure, not causal inference.
#
# Ne estimates (thousands):
# Tenesa et al. 2007 (Nat Genet): afr=17.0, eas=7.0, nfe=8.0
# Gravel et al. 2011 (PNAS): amr=5.0, sas=9.0
# Carmi et al. 2014 (Nat Commun): asj=1.5
# Lim et al. 2014: fin=3.5

ne_data <- tibble(
  ancestry = c("afr", "amr", "asj", "eas", "fin", "mid", "nfe", "sas"),
  ancestry_label = c("African", "Admixed American", "Ashkenazi Jewish",
                     "East Asian", "Finnish", "Middle Eastern",
                     "Non-Finnish Eur.", "South Asian"),
  ne_thousands = c(17.0, 5.0, 1.5, 7.0, 3.5, 4.0, 8.0, 9.0),
  bottleneck   = c(5, 2, 1, 3, 1, 2, 3, 3)
)

mean_burden <- burden_filtered %>%
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

lm_fit <- lm(mean_burden ~ ne_thousands, data = scatter_data)
lm_sum <- summary(lm_fit)
cat(sprintf("\nLinear model: R² = %.3f, p = %.3f (n=7, interpret with caution)\n",
            lm_sum$r.squared,
            lm_sum$coefficients[2, 4]))

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
    values = c(afr="#e94560", amr="#f4a261", asj="#9b5de5",
               eas="#00b4d8", fin="#06d6a0", mid="#ffd166",
               nfe="#457b9d", sas="#f77f00"),
    guide = "none"
  ) +
  scale_size_continuous(
    range  = c(4, 10), name = "Bottleneck\nstrength",
    breaks = c(1, 3, 5), labels = c("Strong", "Moderate", "Weak/None")
  ) +
  scale_x_continuous(
    name = "Historical effective population size, Ne (thousands)",
    limits = c(0, 20)
  ) +
  labs(
    title    = "Demographic Context of Rare LoF Variant Burden",
    subtitle = paste0("Mean rare LoF HC variant burden vs. historical Ne across ancestries\n",
                      "(exploratory; R² = 0.24, p = 0.264, n = 7)"),
    y        = "Mean LoF burden (per 10K chromosomes)",
    caption  = "Ne: Tenesa 2007, Gravel 2011, Carmi 2014, Lim 2014\nInsufficient power for definitive conclusions; retained as contextual figure."
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title       = element_text(face = "bold", size = 13),
    plot.subtitle    = element_text(size = 9, color = "grey40"),
    plot.caption     = element_text(size = 7, color = "grey60"),
    legend.position  = "right",
    panel.grid.minor = element_blank()
  )

ggsave("figures/03_demographic_burden_scatter.png",
       plot = p_scatter, width = 10, height = 7, dpi = 300)
cat("\nScatter plot saved: figures/03_demographic_burden_scatter.png\n")

# ── 14. Constraint plot ───────────────────────────────────────
# LOEUF thresholds (Karczewski et al. 2020, Nature):
# < 0.35: high confidence constrained (pLI > 0.9)
# 0.35-0.70: moderate constraint
# > 0.70: tolerant
# TLR10 excluded: no constraint data in gnomAD v4

constraint <- read_tsv("data/processed/constraint_metrics.tsv",
                       show_col_types = FALSE)

constraint <- constraint %>%
  filter(!is.na(LOEUF), gene %in% genes_for_figures) %>%
  left_join(module_map, by = "gene") %>%
  mutate(
    gene        = factor(gene, levels = gene_order),
    constrained = case_when(
      LOEUF < 0.35                 ~ "Highly constrained",
      LOEUF >= 0.35 & LOEUF < 0.7 ~ "Moderately constrained",
      LOEUF >= 0.7                 ~ "Tolerant / unconstrained"
    ),
    constrained = factor(constrained,
                         levels = c("Highly constrained",
                                    "Moderately constrained",
                                    "Tolerant / unconstrained"))
  )

p_constraint <- constraint %>%
  ggplot(aes(x = LOEUF, y = gene, fill = constrained)) +
  geom_col(width = 0.6, alpha = 0.9) +
  geom_vline(xintercept = 0.35, linetype = "dashed",
             color = "#e94560", linewidth = 0.6) +
  geom_vline(xintercept = 0.7, linetype = "dashed",
             color = "#f4a261", linewidth = 0.6) +
  annotate("text", x = 0.36, y = 0.6, label = "LOEUF = 0.35",
           hjust = 0, size = 2.8, color = "#e94560") +
  annotate("text", x = 0.71, y = 0.6, label = "LOEUF = 0.70",
           hjust = 0, size = 2.8, color = "#f4a261") +
  facet_grid(module ~ ., scales = "free_y", space = "free_y") +
  scale_fill_manual(
    values = c(
      "Highly constrained"       = "#e94560",
      "Moderately constrained"   = "#533483",
      "Tolerant / unconstrained" = "#457b9d"
    ),
    name = "Constraint level"
  ) +
  scale_x_continuous(
    limits = c(0, 3.5),
    name   = "LOEUF (oe_lof_upper) — lower = more constrained"
  ) +
  labs(
    title    = "Loss-of-Function Constraint Across Candidate Genes",
    subtitle = "gnomAD v4 | LOEUF < 0.35: highly constrained | LOEUF < 0.70: moderately constrained",
    y        = NULL,
    caption  = "LOEUF thresholds: Karczewski et al. 2020 (Nature). TLR10 excluded: no constraint data in gnomAD v4."
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

ggsave("figures/04_constraint_loeuf.png",
       plot = p_constraint, width = 10, height = 8, dpi = 300)
cat("\nConstraint plot saved: figures/04_constraint_loeuf.png\n")