# 01_fetch_gnomad.py
# Fetches variant-level data from gnomAD v4 for candidate genes
# using the public GraphQL API. Outputs one TSV per gene.

import requests
import time
import glob
import pandas as pd
from pathlib import Path

GNOMAD_API = "https://gnomad.broadinstitute.org/api"
DATASET    = "gnomad_r4"
OUTPUT_DIR = Path("data/raw")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

GENES = [
    "IKBKB", "IFIH1", "TLR1", "TLR6", "TLR10", "ACKR1", "CD36",
    "STAT4", "IRF5", "TNFAIP3",
    "RUNX2", "COL1A1", "SP7", "LRP5", "VDR",
]

QUERY = """
query GeneVariants($gene_symbol: String!, $dataset: DatasetId!) {
  gene(gene_symbol: $gene_symbol, reference_genome: GRCh38) {
    gene_id
    symbol
    variants(dataset: $dataset) {
      variant_id
      pos
      consequence
      lof
      lof_filter
      lof_flags
      exome {
        ac
        an
        af
        populations {
          id
          ac
          an
        }
      }
      genome {
        ac
        an
        af
        populations {
          id
          ac
          an
        }
      }
    }
  }
}
"""

def fetch_gene(gene_symbol):
    response = requests.post(
        GNOMAD_API,
        json={"query": QUERY, "variables": {"gene_symbol": gene_symbol, "dataset": DATASET}},
        headers={"Content-Type": "application/json"},
        timeout=120,
    )
    response.raise_for_status()
    data = response.json()

    if "errors" in data:
        print(f"  [!] GraphQL errors for {gene_symbol}: {data['errors']}")
        return []

    variants = data.get("data", {}).get("gene", {}).get("variants", [])
    print(f"  -> {len(variants)} variants retrieved")
    return variants


def parse_variants(gene_symbol, variants):
    rows = []
    for v in variants:
        base = {
            "gene":        gene_symbol,
            "variant_id":  v["variant_id"],
            "pos":         v["pos"],
            "consequence": v.get("consequence", ""),
            "lof":         v.get("lof", ""),
            "lof_filter":  v.get("lof_filter", ""),
        }

        for dtype in ("exome", "genome"):
            block = v.get(dtype)
            if not block:
                continue

            rows.append({
                **base,
                "data_type": dtype,
                "ancestry":  "ALL",
                "ac":        block.get("ac", 0),
                "an":        block.get("an", 0),
                "af":        block.get("af", 0.0),
            })

            for pop in block.get("populations", []):
                an = pop.get("an", 0)
                ac = pop.get("ac", 0)
                rows.append({
                    **base,
                    "data_type": dtype,
                    "ancestry":  pop["id"],
                    "ac":        ac,
                    "an":        an,
                    "af":        ac / an if an > 0 else 0.0,
                })

    return pd.DataFrame(rows)


if __name__ == "__main__":
    print(f"Fetching {len(GENES)} genes from gnomAD {DATASET}...\n")

    all_dfs = []
    for gene in GENES:
        print(f"[{gene}]")
        variants = fetch_gene(gene)
        if not variants:
            continue

        df = parse_variants(gene, variants)
        out_path = OUTPUT_DIR / f"{gene}_gnomad_v4.tsv"
        df.to_csv(out_path, sep="\t", index=False)
        print(f"  -> Saved: {out_path}  ({len(df)} rows)\n")

        all_dfs.append(df)
        time.sleep(1)

    # Combine all TSVs including previously downloaded genes
    existing = [pd.read_csv(f, sep="\t") for f in glob.glob("data/raw/*_gnomad_v4.tsv")]
    combined = pd.concat(existing, ignore_index=True)
    combined_path = OUTPUT_DIR / "all_genes_gnomad_v4.tsv"
    combined.to_csv(combined_path, sep="\t", index=False)
    print(f"\nDone. Combined file: {combined_path}")
    print(f"Total rows: {len(combined):,}")