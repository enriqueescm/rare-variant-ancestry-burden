# 03_fetch_constraint.py
# Fetches gnomAD constraint metrics for 15 candidate genes
# LOEUF (oe_lof_upper): lower = more constrained
# pLI: probability of being loss-of-function intolerant (>0.9 = constrained)

import requests
import time
import pandas as pd
from pathlib import Path

GNOMAD_API = "https://gnomad.broadinstitute.org/api"
OUTPUT_DIR = Path("data/processed")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

GENES = [
    "IKBKB", "IFIH1", "TLR1", "TLR6", "TLR10", "ACKR1", "CD36",
    "STAT4", "IRF5", "TNFAIP3",
    "RUNX2", "COL1A1", "SP7", "LRP5", "VDR",
]

QUERY = """
query GeneConstraint($gene_symbol: String!) {
  gene(gene_symbol: $gene_symbol, reference_genome: GRCh38) {
    symbol
    gnomad_constraint {
      pLI
      oe_lof_upper
      oe_mis
      oe_mis_upper
      oe_syn
    }
  }
}
"""

def fetch_constraint(gene_symbol):
    response = requests.post(
        GNOMAD_API,
        json={"query": QUERY, "variables": {"gene_symbol": gene_symbol}},
        headers={"Content-Type": "application/json"},
        timeout=30,
    )
    response.raise_for_status()
    data = response.json()

    if "errors" in data:
        print(f"  [!] Error for {gene_symbol}: {data['errors']}")
        return None

    constraint = data.get("data", {}).get("gene", {}).get("gnomad_constraint", {})
    return constraint


if __name__ == "__main__":
    # Load already fetched genes to avoid repeating API calls
    out_path = OUTPUT_DIR / "constraint_metrics.tsv"
    if out_path.exists():
        existing = pd.read_csv(out_path, sep="\t")
        already_done = set(existing["gene"].tolist())
    else:
        existing = pd.DataFrame()
        already_done = set()

    genes_to_fetch = [g for g in GENES if g not in already_done]
    print(f"Fetching constraint metrics for {len(genes_to_fetch)} genes "
          f"({len(already_done)} already cached)...\n")

    rows = []
    for gene in genes_to_fetch:
        print(f"[{gene}]")
        c = fetch_constraint(gene)
        if not c:
            continue

        pli = c.get("pLI")
        loeuf = c.get("oe_lof_upper")
        pli_str = f"{pli:.3f}" if pli is not None else "NA"
        loeuf_str = f"{loeuf:.3f}" if loeuf is not None else "NA"
        print(f"  pLI={pli_str}  LOEUF={loeuf_str}")

        rows.append({
            "gene":         gene,
            "pLI":          pli,
            "LOEUF":        loeuf,
            "oe_mis":       c.get("oe_mis"),
            "oe_mis_upper": c.get("oe_mis_upper"),
            "oe_syn":       c.get("oe_syn"),
        })
        time.sleep(3)

    # Combine with existing and save
    df_new = pd.DataFrame(rows)
    df_all = pd.concat([existing, df_new], ignore_index=True)
    df_all = df_all[df_all["gene"].isin(GENES)]  # keep only our genes
    df_all.to_csv(out_path, sep="\t", index=False)
    print(f"\nDone. Saved: {out_path}")
    print(df_all.to_string(index=False))

 # ── 14. Constraint plot ───────────────────────────────────────
constraint <- read_tsv("data/processed/constraint_metrics.tsv",
                       show_col_types = FALSE)

constraint <- constraint %>%
  left_join(module_map, by = "gene") %>%
  mutate(
    gene        = factor(gene, levels = gene_order),
    constrained = case_when(
      LOEUF < 0.35                 ~ "Highly constrained",
      LOEUF >= 0.35 & LOEUF < 0.7 ~ "Moderately constrained",
      LOEUF >= 0.7 | is.na(LOEUF) ~ "Tolerant / unconstrained"
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
  annotate("text", x = 0.36, y = 0.6, label = "pLI threshold",
           hjust = 0, size = 2.8, color = "#e94560") +
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
    subtitle = "gnomAD v4 | LOEUF < 0.35: highly constrained | NA = no constraint data",
    y        = NULL,
    caption  = "LOEUF: upper bound of observed/expected LoF ratio confidence interval"
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
  "figures/04_constraint_loeuf.png",
  plot   = p_constraint,
  width  = 10,
  height = 8,
  dpi    = 300
)

cat("\nConstraint plot saved: figures/04_constraint_loeuf.png\n")