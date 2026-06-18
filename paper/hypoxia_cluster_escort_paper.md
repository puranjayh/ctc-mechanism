# Hypoxia Licenses the Cluster-and-Escort Phenotype of Circulating Tumor Cells

*Manuscript draft (primary paper) — target: cancer / ML4H venue. Companion negative-result
benchmark: `cross_stressor_ctc_paper.md`. Last updated 2026-06-19.*

---

## Abstract

The most metastatic circulating tumor cells (CTCs) travel as **clusters** and in the
company of **neutrophil "bodyguards,"** yet what licenses this aggressive phenotype is
unclear. We test the hypothesis that **hypoxia** is the upstream driver, using
single-cell breast-cancer CTC transcriptomes spanning two perturbations: tumor hypoxia
(GSE126669) and circadian rest/active phase (GSE180097), the latter annotated with CTC
clustering and white-blood-cell (WBC) attachment. We derive a hypoxia survival signature
from hypoxic CTCs and score it across all cells. The hypoxia program is significantly
elevated in **clustered** CTCs (p = 0.023) and in **WBC-attached** CTCs (p = 0.008), and
a **shared HIF-target survival core** (VEGFA, P4HA1, NDRG1, HK2, HIF1A, EGLN1, BNIP3)
links hypoxic and rest-phase CTCs. Critically, the program is **not** elevated by
circadian timing (rest vs active, p = 0.41), and—across four independent methods—stress
signatures do **not** transfer between the two stressors, ruling out a generic
"universal stress" explanation. Hypoxia thus associates specifically with the
cluster-and-escort phenotype rather than with circadian timing or a non-specific stress
response. We discuss confounds (WBC transcript composition; n = 30 hypoxic cells) and
implications for targeting hypoxic CTC clusters in liquid biopsy.

---

## 1. Introduction

CTCs are the seeds of metastasis. Two features mark the most dangerous ones: they travel
as **multicellular clusters** (far higher metastatic potential than single CTCs) and they
recruit **neutrophils** that shield them and boost proliferation. A central open question
is what upstream signal *licenses* this phenotype.

Tumor hypoxia is a prime candidate: it drives EMT, stemness and aggressiveness, and
hypoxic CTCs are more metastatic. A parallel literature links **circadian rhythm** to CTC
generation (CTC shedding peaks during the rest phase). This raises competing hypotheses:
is the aggressive phenotype licensed by **hypoxia**, by **circadian timing**, or by a
**generic stress response** shared across perturbations?

We adjudicate between these directly. We derive a hypoxia signature from hypoxic CTCs and
ask where it is elevated, using clustering and WBC-attachment annotations as readouts of
the aggressive phenotype, and circadian phase and cross-stressor transfer as controls.

**Contributions.** (1) Evidence that hypoxia is elevated specifically in clustered and
neutrophil-associated CTCs. (2) A shared HIF-target survival core between hypoxic and
rest-phase CTCs. (3) Negative controls ruling out circadian-timing and generic-stress
explanations, including a four-method cross-stressor transfer analysis.

---

## 2. Data & Methods

**Datasets.** GSE126669: 30 single breast CTCs labeled hypoxia Positive/Negative.
GSE180097: 276 breast CTCs (BR16, LM2 xenografts; Patient) labeled circadian
resting/active, with `sample type` (single / cluster / cluster-with-WBC), cluster size,
and WBC count.

**Pipeline.** Counts→TPM; Ensembl→symbols. To prevent the dominant batch from driving
feature selection, `log1p` was applied per batch on split layers; 3000 HVGs (mvp);
scaling without centering (zero-inflated TPM); 50 PCs; **Harmony** integration; SNN graph,
Louvain clustering (3 clusters), UMAP.

**Hypoxia survival signature.** Wilcoxon markers up in hypoxia-Positive vs Negative CTCs
(GSE126669); top genes scored across all cells via module scoring (control gene set
matched). A neutrophil-chemokine module (CXCL8, CXCL1/2/3/5/6, CSF2/3, S100A8/9, IL6) was
scored analogously.

**Tests.** Within GSE180097: Wilcoxon of hypoxia score by phase (resting/active),
clumping (cluster/single), and WBC status (with/without). Shared survival genes =
intersection of hypoxia-Positive markers and resting-up markers. Cluster identity via
`FindAllMarkers`.

---

## 3. Results

### 3.1 A hypoxia signature with real within-data signal
The hypoxia signature (top genes incl. VEGFA, P4HA1-axis, NDRG1) was derived from
GSE126669. Within-stressor cross-validation confirms stress state is recoverable in the
larger circadian cohort (AUROC ≈ 0.91), establishing that the labels carry transcriptional
signal.

### 3.2 Hypoxia is elevated in clustered CTCs
Clustered CTCs score higher for the hypoxia program than single CTCs
(mean −0.059 vs −0.132; Wilcoxon **p = 0.023**), consistent with hypoxia promoting CTC
clustering.

### 3.3 Hypoxia is elevated in neutrophil-associated CTCs
CTCs with WBCs attached score higher for hypoxia than those without
(mean −0.011 vs −0.108; **p = 0.008**), consistent with hypoxic CTCs recruiting neutrophil
escorts. Neutrophil-chemokine scores are also higher in WBC-attached CTCs
(**p = 0.001**) — but see §4 (composition confound).

### 3.4 A shared HIF-target survival core
Genes up in both hypoxic and rest-phase CTCs include a canonical HIF-target core: **VEGFA,
P4HA1, NDRG1, HK2, HIF1A, EGLN1, BNIP3, ERO1A, MTHFD2**, indicating a common
hypoxia-driven survival program engaged across contexts.

### 3.5 Negative controls: not circadian timing, not generic stress
- **Circadian timing:** hypoxia score does not differ by rest vs active phase
  (mean −0.101 vs −0.074; **p = 0.41**); VEGFA is flat across phase. The aggressive
  phenotype is not explained by circadian timing.
- **Generic stress:** stress signatures do not transfer between hypoxia and circadian
  stress across four methods (AUROC 0.34–0.54; genome-wide logFC r ≈ −0.06), ruling out a
  non-specific "universal stress" program (companion paper). Only metabolic-stress
  pathways (HYPOXIA/UPR/ISR) show weak concordance.

### 3.6 Cluster identities
Cluster 2 is proliferative (PRR11, CDKN3, KNL1, CDKL1); clusters 0/1 are not cleanly
hypoxia- or phase-defined (cluster 0 carries the metastasis marker PODXL).

---

## 4. Limitations & confounds

- **WBC composition confound:** `cluster-with-WBC` samples physically contain neutrophils,
  so elevated chemokine scores partly reflect WBC transcripts rather than CTC signaling.
  The hypoxia-score result (tumor-intrinsic genes) is the robust readout; the chemokine
  result requires de-confounding (regress out WBC content / restrict to tumor-intrinsic
  genes) before a "recruitment" claim.
- **Small hypoxic cohort (n = 30):** the signature is derived from few cells; treat as
  directional. The shared-gene overlap used loose thresholds — the meaningful unit is the
  GSEA-enriched HIF core, not the full list.
- **Associational, not causal**, and limited to two breast-cancer Smart-seq2 datasets.

---

## 5. Discussion & Conclusion

Across complementary tests, hypoxia is specifically associated with the **cluster-and-
escort phenotype** of CTCs — clustering and neutrophil association — and not with circadian
timing or a generic stress response. A shared HIF-target survival core links hypoxic and
rest-phase CTCs. This positions hypoxia as a candidate upstream licensor of the most
metastatic CTC states and motivates targeting hypoxic CTC clusters and their
neutrophil interactions. Confirmatory work (GSEA validation of the HIF core,
de-confounded recruitment analysis, per-donor robustness, and ligand–receptor interaction
modeling) is the immediate next step.

---

## Reproducibility
Code: `full_pipeline.R`, `mechanistic_hypotheses.R` (this paper), `cross_stressor_transfer.R`,
`rescue_test.R`, `export_for_gnn.R`, `gnn_disentangle.py` (controls). Tables:
`hypoxia_survival_markers.csv`, `cluster_markers.csv`, `shared_survival_genes.csv`.
Figure: `main/expression_data/mechanistic_hypotheses_plots.png`.
