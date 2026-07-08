# Ancestry-Stratified Rare LoF Variant Burden in Immune and Skeletal Genes

**Enrique Estévez Campo** · Biological Anthropologist & PhD in Biomedicine  
---

## Overview

This project investigates how demographic history shapes the landscape of rare
loss-of-function (LoF) genetic variation across human populations, using a
curated panel of 15 candidate genes spanning three biological modules: innate
immunity, autoimmunity, and skeletal biology.

The central hypothesis is that the same evolutionary forces that diversified
human morphology — genetic drift, founder effects, and pathogen-driven
selection — leave distinct signatures in the spectrum of rare functional
variants across ancestries. This analysis bridges biological anthropology and
population genomics, connecting macroscopic patterns of human diversity to
their molecular underpinnings.

Inspired by recent work at deCODE Genetics identifying ancestry-specific rare
variants with large effects on disease risk (Thorlacius et al. 2025), this
project asks: **what does the landscape of rare LoF variation look like when
ancestry is treated as a primary axis of analysis?**

---

## Biological modules

| Module | Genes | Rationale |
|--------|-------|-----------|
| A: Innate immunity | IKBKB, IFIH1, TLR1, TLR6, TLR10, ACKR1, CD36 | Genes under differential positive selection across populations; known ancestry-specific variants |
| B: Autoimmunity | STAT4, IRF5, TNFAIP3 | Trans-ancestral GWAS loci for lupus and autoimmune disease; interferon axis |
| C: Skeletal biology | RUNX2, COL1A1, SP7, LRP5, VDR | Osteogenesis and bone density genes; connect to author's background in skeletal biology |

---

## Data

- **Source**: gnomAD v4.1 (Karczewski et al. 2020; Collins et al. 2024)
- **Dataset**: Exomes only (734,947 individuals)
- **Ancestry groups**: African, Admixed American, Ashkenazi Jewish, East Asian,
  Finnish, Middle Eastern, Non-Finnish European, South Asian
- **Variant filter**: MAF < 1% per ancestry | LoF high-confidence (LOFTEE HC)
- **Burden metric**: Σ(AC) / max(AN) × 10,000 per ancestry per gene

All data retrieved via the gnomAD GraphQL public API. No individual-level data
used. Full pipeline is reproducible from scripts in this repository.

---

## Key findings

**CD36 shows the strongest ancestry-differentiated LoF burden**, with African
(456 per 10K chr) and East Asian (422) populations carrying substantially
higher burden than Finnish (6) or Non-Finnish European (52). This is consistent
with documented positive selection at CD36 driven by differential malaria
exposure and is one of the best-characterized examples of pathogen-driven
adaptation in the human genome.

**IFIH1 shows broad elevation across multiple non-European ancestries**,
including Admixed American (397), Middle Eastern (253), and Non-Finnish
European (350). IFIH1 encodes MDA5, a cytosolic RNA sensor involved in antiviral
innate immunity, and has been previously identified as a target of positive
selection in African populations.

**Module B (autoimmunity) genes show uniformly low LoF burden** across all
ancestries, consistent with strong purifying selection. TNFAIP3, a negative
regulator of NF-κB, has the lowest LOEUF (0.217) of all candidate genes.

**Module C (skeletal) genes are highly constrained**: COL1A1 (LOEUF = 0.155)
and RUNX2 (LOEUF = 0.328) are among the most intolerant to LoF variation in
the human genome, consistent with their essential roles in osteogenesis.

**The demographic scatter (Ne vs. burden) shows a positive trend but does not
reach statistical significance** (R² = 0.24, p = 0.264, n = 7), likely due
to insufficient power with only 7 ancestry groups and gene-specific variation
in LoF constraint. This is documented transparently as a limitation.

---

## Figures

| Figure | Description |
|--------|-------------|
| `figures/01_burden_heatmap.png` | LoF burden heatmap across 14 genes × 8 ancestries |
| `figures/02_ancestry_specific_lollipop.png` | Normalized ancestry-specific LoF variants per gene |
| `figures/03_demographic_burden_scatter.png` | Demographic context: Ne vs. mean LoF burden (exploratory) |
| `figures/04_constraint_loeuf.png` | LOEUF constraint across candidate genes |

---

## Repository structure

rare-variant-ancestry-burden/
├── data/
│   ├── raw/                    # Per-gene gnomAD v4 TSVs (fetched via API)
│   └── processed/              # Burden tables, constraint metrics
├── scripts/
│   ├── 01_fetch_gnomad.py      # gnomAD GraphQL API data retrieval
│   ├── 02_burden_analysis.R    # Filtering, burden calculation, figures
│   └── 03_fetch_constraint.py  # LOEUF/pLI constraint metrics
├── figures/                    # Publication-quality PNG figures
└── README.md

---

## Methodological notes

- **LoF HC only**: missense variants excluded from primary burden analysis.
  Without CADD pathogenicity filtering, missense variants inflate burden with
  likely benign variants (Pearson r = 0.46 between missense+LoF and LoF-only
  burden, confirming poor concordance). This is a deliberate conservative choice.

- **Burden denominator**: max(AN) per gene × ancestry used as denominator,
  approximating total callable chromosomes and avoiding downward bias from
  low-coverage variants.

- **MAF filter**: applied per-ancestry (not globally), capturing variants rare
  within each specific population — the biologically relevant unit for
  ancestry-stratified analysis.

- **TLR10**: excluded from burden figures — no LoF HC variants in gnomAD v4
  exomes for this gene.

- **ACKR1**: retained in analysis; zero LoF HC burden across all ancestries,
  consistent with its high LOEUF (3.06) indicating tolerance to LoF variation.

---

## Portfolio context

This is the third project in a portfolio connecting biological anthropology to
population genomics:

1. [morphological-echoes-genome](https://github.com/enriqueescm/morphological-echoes-genome) —
   Population structure and selection signals from 1000 Genomes Phase 3
2. [ilium-trabecular-microstructure](https://github.com/enriqueescm/ilium-trabecular-microstructure) —
   Trabecular bone microstructure from original microCT data (Granada Osteological Collection)
3. **This project** — Ancestry-stratified rare LoF variant burden in immune and skeletal genes

The three projects share a common thread: the same demographic forces that
shaped observable human morphological diversity also structure genetic variation
at the population level — from gross anatomy, to bone microstructure, to the
spectrum of rare functional variants.

---

## References

- Karczewski et al. (2020). The mutational constraint spectrum quantified from
  variation in 141,456 humans. *Nature*, 581, 434–443.
- Collins et al. (2024). A cross-disorder dosage sensitivity map of the human
  genome. *Cell*, 187, 1–15. [gnomAD v4]
- Tenesa et al. (2007). Recent human effective population size estimated from
  linkage disequilibrium. *Genome Research*, 17, 520–526.
- Gravel et al. (2011). Demographic history and rare allele sharing among human
  populations. *PNAS*, 108, 11983–11988.
- Carmi et al. (2014). Sequencing an Ashkenazi reference panel supports
  population-targeted personal genomics. *Nature Communications*, 5, 4835.
- Lim et al. (2014). Distribution and medical impact of loss-of-function
  variants in the Finnish founder population. *PLOS Genetics*, 10, e1004494.
- Thorlacius et al. (2025). African-ancestry-specific variant IKKβ p.Glu502Lys
  confers high lupus risk. *Nature Genetics*, 57, 2980–2986.

---

## Author

**Enrique (Kike) Estévez Campo**  
PhD in Biomedicine · Biological Anthropologist  
Population genomics | Human genetic variation
[GitHub](https://github.com/enriqueescm) · [LinkedIn](https://www.linkedin.com/in/enrique-estevez-campo/)