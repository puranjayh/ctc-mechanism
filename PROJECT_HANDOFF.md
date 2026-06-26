# CTC Mechanism Project — Session Handoff

> ⚠️ **Partly superseded.** §6 below (mechanistic hypotheses) reported results later shown
> to be **donor/contamination artifacts**. See `README.md` and `deconfound_results.md` for
> the corrected findings. This file is kept for background/provenance.
>
> Paste this whole file into a new Claude session on the new laptop to continue with full context.
> Last updated: 2026-06-19. Repo: https://github.com/puranjayh/ctc-mechanism

---

## 1. What this project is

Single-cell RNA-seq study of **circulating tumor cells (CTCs)** in breast cancer,
integrating two datasets to probe how CTCs survive in the bloodstream. The work
has evolved through three framings (in order):

1. **Original idea:** find a common "stress" gene program shared across stressors
   → a universal CTC stress signature / liquid-biopsy biomarker.
2. **ML framing (for NeurIPS workshops):** test *cross-stressor transfer* — does a
   stress signature learned on one stressor generalize to another? (Result: NO —
   see §5.) Reframed as a benchmark + negative-result paper.
3. **Current framing (PhD's direction):** a **mechanistic model of metastasis** —
   hypoxia DRIVES CTCs to (a) cluster and (b) recruit neutrophil "bodyguards."
   This is the active direction. (Results: §6.)

---

## 2. Datasets

| Domain | Accession | Stressor / variable | Cells | Notes |
|---|---|---|---|---|
| A | GSE126669 | Hypoxia (Positive/Negative) | 30 | single CTCs, breast; models Br61/Br16 |
| B | GSE180097 | Circadian (resting/active) | 276 (of 306) | BR16/LM2 xenograft + Patient; Aceto-style circadian CTC data |

**GSE180097 metadata fields** (`main/metadata/metadata_GSE180097.csv`):
`donor` (BR16 138 / LM2 132 / Patient 36), `timepoint` (resting 143 / active 139 +
24 numeric rows to drop), `sample type` (ctc_cluster 124 / ctc_single 117 /
**ctc_cluster_wbc 65** ← these have WBCs attached), `number of ctcs` (cluster size),
`number of wbcs attached` (65 non-null).
**GSE126669 metadata** (`main/metadata/metadata_GSE126669.csv`): `hypoxia` Pos/Neg.

**Unified label:** `stress_condition` = Stressed/Baseline, where resting=Stressed,
active=Baseline, hypoxia Positive=Stressed, Negative=Baseline.

---

## 3. Pipeline (decided with PhD)

- Counts → TPM; Ensembl → gene symbols (org.Hs.eg.db).
- **log1p applied PER BATCH on split layers** (NOT after JoinLayers) — so HVG
  selection isn't dominated by one batch. (PhD insisted on this.)
- HVGs via **mvp** method, 3000 genes (safer for Smart-seq2 TPM).
- **ScaleData(do.center = FALSE)** — protects zero-inflated TPM.
- RunPCA(npcs = 50) → **Harmony** integration (batch var = `dataset_id`).
- FindNeighbors/FindClusters(res 0.5) → 3 clusters; RunUMAP (`umap.harmony`).
- **CCA integration was abandoned** (PhD said harmony only).

---

## 4. Key files in repo

**Pipeline / analysis (R):**
- `full_pipeline.R` — full preprocessing + Harmony (harmony-only version).
- `cross_stressor_transfer.R` — gene-level transfer test (Method 1).
- `rescue_test.R` — logFC concordance + pathway transfer (Methods 2–3).
- `export_for_gnn.R` — exports node_features.csv / edges.csv / node_meta.csv.
- `mechanistic_hypotheses.R` — **current direction**: H1/H2 tests (§6).

**GNN (Python, pure PyTorch, CPU):**
- `gnn_disentangle.py` — disentanglement GCN (Method 4).

**Paper:**
- `paper/hypoxia_cluster_escort_paper.md` + `.docx` — **PRIMARY manuscript** (mechanistic:
  hypoxia → CTC clustering + neutrophil escort; includes the negative results as controls).
- `paper/cross_stressor_ctc_paper.md` + `.docx` — companion workshop manuscript
  (cross-stressor negative-result benchmark + GNN).
- `paper/build_docx_mechanistic.js`, `paper/build_docx.js` — regenerate the Word docs.

**Result tables (CSV):**
- `hypoxia_survival_markers.csv` — genes up in hypoxia-Positive CTCs.
- `cluster_markers.csv` — FindAllMarkers per cluster.
- `shared_survival_genes.csv` — genes shared by hypoxic AND resting CTCs (847; loose).
- `concordant_stress_signature.csv` — per-gene logFC in both stressors.

**NOT in git (gitignored, too big — transfer via cloud/USB):**
- All `*.RData` (180–775 MB), `*.rds`, `*.png`, and the `expression_data/` raw data.
- Latest workspace: `post_mechanistic_hypotheses.RData` (run `save.image(...)` to create).
- NOTE: `main/` is a **nested git repo** — files inside it (incl. metadata, figures)
  are tracked separately, not by the parent repo.

---

## 5. Cross-stressor transfer — NEGATIVE (settled, don't re-run)

Question: does a stress classifier trained on one stressor predict the other?
Four independent methods all say **no shared transferable program**:

| Method | Hyp→Circ | Circ→Hyp | Note |
|---|---|---|---|
| L1-logistic (genes) | 0.504 | 0.540 | perm p = 0.41 / 0.68; 0 shared genes |
| Pathway-score transfer | 0.488 | 0.446 | — |
| Disentanglement GNN | 0.492±0.08 | 0.342±0.11 | below chance one way |
| Genome-wide logFC corr | r=−0.059 | ρ=−0.070 | negligible effect |

Within-domain signal exists (circadian CV AUROC ≈ 0.91). Only positive nuance:
**pathway-level concordance limited to metabolic stress** (HYPOXIA/UPR/ISR move same
direction in both; NF-κB/p53/dormancy are stressor-specific).
**Conclusion:** CTC stress is stressor-specific, not universal. (Caveat: hypoxia n=30
is underpowered.) **Do not keep tuning the GNN — that would be p-hacking.**

---

## 6. Mechanistic hypotheses — CURRENT, partially POSITIVE

PhD's model: hypoxia → CTCs cluster → recruit neutrophil bodyguards.
(`mechanistic_hypotheses.R`; hypoxia signature learned from GSE126669, scored on all cells.)

| Hypothesis | Result | Verdict |
|---|---|---|
| H1 timing: hypoxia ↑ at rest | resting −0.10 vs active −0.07, **p=0.41** | ❌ NOT supported |
| H1b: hypoxia ↑ in clusters | cluster −0.06 vs single −0.13, **p=0.023** | ✅ supported |
| H2: hypoxia ↑ with WBCs attached | with_WBC −0.01 vs no_WBC −0.11, **p=0.008** | ✅ supported |
| H2b: neutrophil chemokines ↑ with WBCs | +1.61 vs −0.31, **p=0.001** | ⚠️ supported but CONFOUNDED |

**Interpretation:** hypoxia tracks the **aggressive cluster + immune-escort phenotype**,
NOT circadian timing. The "night-time" part of the model is not supported here.
**Caveats:** (1) H2b chemokines confounded — `ctc_cluster_wbc` samples physically
contain neutrophils, so CXCL8/S100A8/9 partly = WBC transcripts, not CTC "calling out";
the hypoxia-score result (H2a) is the trustworthy one. (2) `shared_survival_genes.csv`
(847 genes) is too loose, but contains a real HIF-target core: VEGFA, P4HA1, NDRG1,
HK2, HIF1A, EGLN1, BNIP3, ERO1A, MTHFD2. (3) Cluster 2 = proliferating (PRR11/CDKN3/
KNL1/CDKL1), not stress-defined.

---

## 7. Environment

- Windows; R 4.6 (RStudio); R libs: Seurat v5, ggplot2, patchwork, org.Hs.eg.db,
  AnnotationDbi, SingleCellExperiment, glmnet, pROC. (Rtools NOT installed — install
  binary packages only.)
- Python 3.12 with torch 2.12, scikit-learn, pandas, numpy (CPU; no torch_geometric needed).
- Git user: Puranjayh; remote: https://github.com/puranjayh/ctc-mechanism (main branch).
  NOTE: remote was once force-rebuilt via web uploads — local `full_pipeline.R` differs
  from the remote copy; we deliberately did not overwrite it.

---

## 8. Resume on the new laptop

```bash
git clone https://github.com/puranjayh/ctc-mechanism
```
Then download from cloud/USB (NOT in git): `ctc_obj.rds` (or
`post_mechanistic_hypotheses.RData`) + the `expression_data/` raw data folder.
In R:
```r
obj <- readRDS("ctc_obj.rds")   # or load("post_mechanistic_hypotheses.RData")
```

---

## 9. Next steps (recommended order)

1. **GSEA-validate the HIF-target survival core** from `shared_survival_genes.csv`
   (use fgsea/Hallmark; expect HALLMARK_HYPOXIA strongly enriched). Report the tight
   core, not all 847 genes.
2. **De-confound H2b:** restrict the neutrophil-chemokine test to genes that are
   tumor-intrinsic, or regress out a WBC-content score, to show CTCs (not the attached
   WBCs) drive recruitment signaling.
3. **Strengthen H1b/H2:** report effect sizes + per-donor breakdown (BR16/LM2/Patient)
   so it's not a single-model artifact.
4. **Decide the paper:** either the mechanistic story ("hypoxia → CTC clustering +
   neutrophil escort") for a cancer/ML4H venue, or the cross-stressor negative-result
   benchmark for a NeurIPS workshop. The mechanistic one is now the stronger, more
   positive result.

---

## 10. One-paragraph context to paste into a new Claude

> I'm continuing a single-cell RNA-seq project on circulating tumor cells (CTCs) in
> breast cancer, integrating GSE126669 (hypoxia, 30 cells) and GSE180097 (circadian
> resting/active, 276 cells; some CTCs have neutrophils/WBCs attached) via Seurat v5 +
> Harmony. We established that cross-stressor stress signatures do NOT transfer (4
> methods, all ~chance — settled negative result). My PhD's current direction is a
> mechanistic model: hypoxia drives CTC clustering and neutrophil "bodyguard"
> recruitment. We found hypoxia score is higher in clustered CTCs (p=0.023) and
> WBC-attached CTCs (p=0.008) but NOT higher at rest (p=0.41), and a shared HIF-target
> survival core (VEGFA/P4HA1/NDRG1/HK2/HIF1A) between hypoxic and rest CTCs. Next I want
> to GSEA-validate that core and de-confound the neutrophil-chemokine signal. Code is in
> the repo (mechanistic_hypotheses.R etc.); the Seurat object is in ctc_obj.rds. See
> PROJECT_HANDOFF.md for full detail.
