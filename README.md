# Hypoxia, circadian timing, and confounding in circulating tumor cells

A single-cell RNA-seq re-analysis of breast-cancer **circulating tumor cells (CTCs)**,
integrating two public datasets to ask what makes CTCs aggressive in the bloodstream —
and, just as importantly, a worked example of how easily **dataset confounding can
manufacture a convincing-but-false mechanism.**

> **TL;DR** — We set out to test a clean hypothesis: that low oxygen (hypoxia) drives CTCs
> to cluster and to recruit neutrophil "bodyguards." Early results looked supportive.
> Proper de-confounding then showed two of the three positives were **artifacts** — of
> cell-line composition (Simpson's paradox) and of immune-cell contamination. What
> survives is honest and clean: a **GSEA-validated hypoxia survival program**, and the
> independent finding that **rest-phase CTCs are more proliferative** (reproducing the
> published "cancer spreads during sleep" biology). The corrected story is stronger —
> and more defensible — than the original.

---

## Background

CTCs are tumor cells that have entered the blood; they are the seeds of metastasis. Two
candidate drivers of CTC aggressiveness are **hypoxia** (low oxygen) and **circadian
timing** (the time of day a cell is shed). We integrated two breast-cancer CTC datasets to
probe both.

| Dataset | Accession | Variable | Cells |
|---|---|---|---|
| Hypoxia | [GSE126669](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE126669) | Hypoxia Positive/Negative | 30 single CTCs |
| Circadian | [GSE180097](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE180097) | Resting/active phase; some CTCs have WBCs attached | 276 CTCs (BR16 / LM2 / Patient) |

Pipeline: counts → TPM → log-normalization → HVG selection → PCA → **Harmony** integration
(by dataset) → clustering/UMAP. Hypoxia and proliferation programs are scored per cell with
`AddModuleScore`; enrichment is tested with `fgsea` against MSigDB Hallmark.

---

## Key results (current, de-confounded)

| Result | Verdict | Evidence |
|---|---|---|
| **Hypoxia survival core is real** | ✅ validated | `fgsea`: `HALLMARK_HYPOXIA` enriched (padj = 0.02); the curated HIF core (VEGFA, P4HA1, NDRG1, HK2, …) is the #1 enriched set |
| **Rest-phase CTCs are more proliferative** | ✅ robust | Donor-controlled p = 0.0036; same direction in all 3 donors. Consistent with the published sleep-accelerates-metastasis result |
| Hypoxia higher in *clustered* CTCs | ❌ artifact | Pooled p = 0.023, but **vanishes within every donor** (Simpson's paradox); donor-adjusted p = 0.56 |
| Hypoxia higher in *WBC-attached* CTCs | ❌ artifact | Pooled p = 0.008 → donor-adjusted p = 0.64 |
| CTCs "recruit" neutrophils (chemokines) | ❌ contamination | The signal is the single gene `S100A9` (a neutrophil-intrinsic protein); recruiter genes (CSF3/CSF2/CXCL1–8) aren't even detected |
| Hypoxia tracks circadian phase | ❌ not supported | p = 0.41 (flat across rest/active) |
| Universal cross-stressor stress signature | ❌ negative | 4 independent methods, all ≈ chance — CTC stress is stressor-specific |

Full numbers and the per-donor breakdown: **[`deconfound_results.md`](deconfound_results.md)**.

**Figures:**
- `deconfound_simpson.png` — the hypoxia–WBC association (pooled "illusion" vs flat within each donor).
- `prolif_timing.png` — proliferation higher at rest in every donor (the robust result).
- `gsea_hypoxia_enrichment.png` — the GSEA enrichment plot for `HALLMARK_HYPOXIA`.

---

## What this project is (and isn't)

It is a careful, honest re-analysis that **reproduces known biology** (the hypoxia program;
sleep-phase proliferation) and **demonstrates how confounders fake mechanisms** in
integrated single-cell data. It is **not** a novel-mechanism discovery — and a key part of
the work was recognizing that and correcting course. The de-confounding is the most useful
contribution here.

---

## Repository structure

**Current analysis (the de-confounding re-analysis):**
| File | What it does |
|---|---|
| `gsea_deconfound.R` | GSEA-validates the hypoxia core; de-confounds the neutrophil signal |
| `perdonor_breakdown.R` | Per-donor stratification — the Simpson's-paradox proof |
| `metastasis_timing.R` | Tests rest-phase proliferation / clustering (donor-controlled) |
| `discover_immune.R`, `inspect_obj.R` | Helper scripts (gene availability, object structure) |
| `*_figure.R` | Generate the result figures |
| `deconfound_results.md` | Detailed write-up of all results |

**Pipeline & setup:**
| File | What it does |
|---|---|
| `full_pipeline.R` | Preprocessing + Harmony integration → builds the Seurat object |
| `install_packages.R` | One-time install of R dependencies |
| `mechanistic_hypotheses.R` | Original hypothesis tests (pre-de-confounding; kept for provenance) |

**Earlier work — cross-stressor transfer (settled negative result):**
| File | What it does |
|---|---|
| `cross_stressor_transfer.R`, `rescue_test.R` | Gene- and pathway-level transfer tests |
| `gnn_disentangle.py`, `export_for_gnn.R` | Disentanglement GNN (PyTorch) + its data export |

**Data & docs:** `metadata/` (GEO metadata + parsers), `expression_data/` (QC figures + loader).

---

## Reproduce

```bash
# 1. install R dependencies (one-time)
Rscript install_packages.R

# 2. run the de-confounding analyses (needs the processed object, see Data availability)
Rscript gsea_deconfound.R
Rscript perdonor_breakdown.R
Rscript metastasis_timing.R
```

Requires R ≥ 4.4 with Seurat v5, fgsea, msigdbr, org.Hs.eg.db. Python parts use PyTorch (CPU).

## Data availability

Raw data are public (GEO accessions above). The **processed Seurat object** (`ctc_obj.rds`,
~53 MB) and large intermediate files are **not** committed (size); rebuild it from the GEO
data via `full_pipeline.R`, or open an issue to request it. Patient samples derive from the
public GSE180097 deposition.

## Limitations

Small cohorts (30 hypoxia cells; 276 circadian CTCs across 3 donors), two breast datasets
only, associational (not causal). Effects are reported with donor controls and per-donor
breakdowns precisely because the data are confounder-prone.

## Acknowledgments

Built on the public GSE126669 and GSE180097 depositions (the latter from the Aceto-lab
sleep/metastasis study). Developed as a mentored research project. See `LICENSE`.
