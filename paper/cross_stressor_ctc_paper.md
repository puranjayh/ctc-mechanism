# Stressor-Specific, Not Universal: A Cross-Stressor Transfer Benchmark for Circulating Tumor Cells

*Workshop manuscript draft — target: NeurIPS LMRL / ML4H (negative-results & benchmarks track)*

---

## Abstract

Circulating tumor cells (CTCs) endure multiple physiological stresses in the
bloodstream, motivating a popular assumption: that CTCs activate a *conserved,
universal stress program* that could serve as a liquid-biopsy biomarker. This
assumption is rarely tested directly. We introduce a **cross-stressor transfer
benchmark** that operationalizes it as a falsifiable generalization question —
*does a stress signature learned from one stressor predict a biologically
orthogonal stressor?* — and a **disentanglement graph neural network (GNN)**
that explicitly separates shared from stressor-specific transcriptional axes.
Applying four independent methods (sparse-linear transfer, genome-wide
fold-change concordance, pathway-activity transfer, and the disentanglement GNN)
to single-cell breast-cancer CTCs under hypoxic stress (GSE126669) and
circadian-rest stress (GSE180097), we find that **cross-stressor transfer does
not exceed chance** (AUROC 0.34–0.54; concordance r ≈ −0.06). The CTC stress
response is largely **stressor-specific, not universal**; the only shared
component is a weak concordance confined to metabolic-stress pathways
(hypoxia/UPR/ISR), while inflammatory, p53, and dormancy programs diverge. A
sensitivity analysis shows the result is robust to, though partially limited by,
the small hypoxic cohort (n = 30). Our benchmark and the convergent negative
finding caution against single-signature stress biomarkers for CTCs and provide
a reusable protocol for testing biological generalization in rare-cell regimes.

---

## 1. Introduction

CTCs are the seeds of metastasis and the substrate of liquid biopsy. In
circulation they face hypoxia, shear stress, anchorage loss, oxidative stress
and immune attack, and are widely hypothesized to mount adaptive stress
responses that promote survival and dormancy. A natural and attractive
corollary — increasingly assumed in the liquid-biopsy literature — is that these
responses converge on a **universal stress signature** detectable across
contexts and usable as a biomarker.

We argue this corollary is a *generalization claim* and should be tested as one.
If a universal stress program exists, a classifier trained to recognize one
stressor should transfer, zero-shot, to a different stressor. We turn two
orthogonal perturbations — **hypoxia** and **circadian rest phase** — into source
and target domains and ask whether stress transfers between them.

**Contributions.**
1. **A cross-stressor transfer benchmark** for CTCs: a leave-one-stressor-out
   protocol with permutation and power controls, evaluated by four complementary
   methods.
2. **A disentanglement GNN** that factors each cell's embedding into a
   domain-invariant *shared-stress* axis and a *stressor-specific* axis, the
   first method (to our knowledge) aimed at isolating a shared stress program in
   rare cells.
3. **A convergent negative result**: cross-stressor transfer fails across all
   four methods, overturning the universal-stress assumption and surfacing a
   weaker, pathway-restricted concordance as a hypothesis for future work.

---

## 2. Related Work

**Universal stress responses.** Conserved multi-stress programs are classic in
biology (e.g., the yeast Environmental Stress Response) and have inspired
literature-derived stress-response gene sets in mammalian toxicology and
multi-stress classifiers in plants (StressGenePred). These establish the
*concept* of shared stress programs but do not test cross-stressor transfer in
CTCs.

**CTC stress and heterogeneity.** Single-cell studies show CTCs activate
stress-tolerance and immune-evasion programs along the vascular journey and that
hypoxia increases CTC aggressiveness and stemness. None test whether distinct
stressors share a *transferable* program.

**Graph and transfer learning for single cells.** GNN-based label transfer
across single-cell datasets (scGCN) and CTC-specific transfer learning
(CTC-Tracer, correcting primary→CTC shift) are established. We repurpose the
graph-transfer machinery for a *generalization-testing* rather than
label-imputation goal, and add explicit shared/specific disentanglement.

---

## 3. Data

| Domain | Accession | Stressor | Cells | Stressed / Baseline |
|---|---|---|---|---|
| A | GSE126669 | Hypoxia (Positive/Negative) | 30 | 14 / 16 |
| B | GSE180097 | Circadian (resting=Stressed / active=Baseline) | 276 | 139 / 137 |

Both are Smart-seq2 breast-cancer CTC datasets (GSE180097 spans BR16, LM2
xenograft and patient samples). We define a unified binary label
`stress ∈ {Stressed, Baseline}` per the biology (resting phase and hypoxia
positive = Stressed).

---

## 4. Methods

**Preprocessing & integration.** Counts were converted to TPM and Ensembl IDs
mapped to symbols. To avoid the dominant batch driving feature selection,
`log1p` was applied *per batch* on split layers; highly variable genes were
selected with the mean-variance-plot method (3000 genes); data were scaled
without centering (appropriate for zero-inflated TPM); 50 PCs were computed and
batches integrated with **Harmony**. A shared-nearest-neighbor (SNN) cell graph
was built on the Harmony embedding.

**Transfer protocol.** For each ordered pair (source, target) ∈
{(A,B),(B,A)} we train on the source stressor's labels and evaluate zero-shot on
the target, reporting AUROC. Controls: (i) within-domain cross-validation to
confirm a within-stressor signal exists; (ii) a **label-permutation null**
(200 shuffles) for significance; (iii) a **power/sensitivity** discussion given
n_A = 30.

**Method 1 — Sparse-linear transfer.** L1-logistic regression on HVG expression.

**Method 2 — Fold-change concordance.** Per-gene Wilcoxon log2FC
(Stressed vs Baseline) computed within each domain; the two genome-wide log2FC
vectors are correlated (Pearson/Spearman).

**Method 3 — Pathway-activity transfer.** Cells scored for six curated stress
programs (HYPOXIA, UPR, P53, NF-κB, ISR, DORMANCY) via module scoring; transfer
is evaluated on the 6-dimensional pathway-score representation, and per-pathway
directional concordance is reported.

**Method 4 — Disentanglement GNN.** A two-layer graph convolutional encoder over
the SNN graph maps each cell to an embedding split into a *shared* and a
*specific* block. Training combines (a) stress prediction on the source domain
from both blocks, (b) a **gradient-reversal domain adversary** forcing the shared
block to be dataset-invariant, and (c) an **orthogonality penalty** decorrelating
shared and specific blocks. Only the shared block is used for cross-stressor
prediction. Reported over 10 seeds (mean ± sd).

---

## 5. Results

### 5.1 A within-stressor signal exists
Within-domain cross-validation recovers stress state in the larger circadian
cohort (CV AUROC ≈ **0.91**), confirming the labels carry real transcriptional
signal. (The hypoxic cohort, n = 30, is too small for a stable within-domain
estimate — see §6.)

### 5.2 Cross-stressor transfer fails across all four methods

| Method | Hypoxia → Circadian | Circadian → Hypoxia | Significance |
|---|---|---|---|
| Sparse-linear (genes) | 0.504 | 0.540 | perm p = 0.41 / 0.68 |
| Pathway-score transfer | 0.488 | 0.446 | — |
| Disentanglement GNN | 0.492 ± 0.075 | 0.342 ± 0.106 | — |
| Genome-wide log2FC concordance | Pearson r = −0.059 | Spearman r = −0.070 | negligible effect (p driven by n=10,564) |

No method exceeds chance; the GNN's Circadian→Hypoxia direction is *below*
chance, echoing the slightly negative genome-wide concordance — the two stress
responses are, if anything, mildly opposed at the single-cell level. The
sparse-linear model selected **zero** genes shared and sign-concordant across
both domains.

### 5.3 The one shared component is pathway-restricted
Per-pathway directional analysis shows the only concordant programs are the
metabolic/energy-stress axis:

| Pathway | Hypoxia shift | Circadian shift | Concordant? |
|---|---|---|---|
| HYPOXIA | +0.519 | +0.136 | ✅ |
| UPR | +0.009 | +0.196 | ✅ |
| ISR | +0.006 | +0.148 | ✅ |
| P53 | −0.120 | +0.148 | ✗ |
| NF-κB | +0.700 | −0.096 | ✗ |
| DORMANCY | +0.220 | −0.038 | ✗ |

A shared adaptive *metabolic*-stress core (hypoxia/UPR/ISR) moves in the same
direction under both stressors, but it is too weak to support transfer, and
inflammatory, p53 and dormancy responses are stressor-specific.

---

## 6. Limitations

The hypoxic cohort is small (n = 30): both training on it (Hypoxia→Circadian) and
predicting it (the unstable 0.342 ± 0.106) are underpowered, so part of the null
reflects statistical power, not only biology. The rigorous claim is therefore
*"with available data, a shared transferable CTC stress axis is not detectable."*
We deliberately did **not** tune the GNN until a metric crossed an arbitrary
threshold, to avoid optimistic bias; the reported configuration is fixed a
priori. Both datasets are breast-cancer Smart-seq2; generality to other cancers
and platforms is untested. Pathway concordance is hypothesis-generating, not
confirmatory.

---

## 7. Conclusion

Across four independent methods, hypoxic and circadian-rest stress responses in
CTCs do not share a transferable signature: stress in CTCs is **stressor-specific,
not universal**, with only a weak, pathway-restricted metabolic-stress overlap.
Beyond the finding, we contribute a reusable cross-stressor transfer benchmark
and a disentanglement GNN for probing biological generalization in rare-cell
regimes. Our results caution against single-signature stress biomarkers for
liquid biopsy and motivate (i) larger, matched multi-stressor CTC cohorts and
(ii) methods that model stressor-specific, rather than universal, adaptation.

---

## Reproducibility

All code is released: `full_pipeline.R` (preprocessing + Harmony),
`cross_stressor_transfer.R` (Method 1), `rescue_test.R` (Methods 2–3),
`export_for_gnn.R` + `gnn_disentangle.py` (Method 4). Result tables:
`concordant_stress_signature.csv`.
