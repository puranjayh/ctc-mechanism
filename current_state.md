# Current State — CTC Mechanism Project

> Snapshot of what's been proven / disproven, what to do next, and the planned paper.
> Last updated: 2026-06-19. Companion to `PROJECT_HANDOFF.md`.

---

## 1. Scoreboard — what's proven, disproven, open

### ✅ Supported by the data
| Finding | Evidence | Strength |
|---|---|---|
| Hypoxia program is higher in **clustered** CTCs (vs single) | cluster −0.06 vs single −0.13, Wilcoxon **p = 0.023** | Moderate |
| Hypoxia program is higher in **WBC-attached** CTCs | with_WBC −0.01 vs no_WBC −0.11, **p = 0.008** | Moderate |
| A **shared HIF-target survival core** exists between hypoxic & rest CTCs | VEGFA, P4HA1, NDRG1, HK2, HIF1A, EGLN1, BNIP3, ERO1A, MTHFD2 present in overlap | Promising, needs GSEA |
| Within-stressor signal is real | circadian CV AUROC ≈ 0.91 | Strong |
| Pathway-level concordance limited to metabolic stress | HYPOXIA/UPR/ISR concordant; NF-κB/p53/dormancy not | Moderate |

### ❌ Disproven / not supported
| Claim | Evidence | Status |
|---|---|---|
| **Universal cross-stressor stress signature** (transfers between stressors) | 4 methods all ≈ chance (AUROC 0.34–0.54; logFC r ≈ −0.06; 0 shared genes) | Settled NO |
| **Night-time timing** — hypoxia higher during rest phase | resting −0.10 vs active −0.07, **p = 0.41**; VEGFA flat by phase | NOT supported |
| Disentanglement GNN recovers a transferable shared axis | 0.49 / 0.34 AUROC, does not beat baseline | NO (do not keep tuning) |

### ⚠️ Supported but CONFOUNDED (cannot claim yet)
| Claim | Issue |
|---|---|
| Hypoxic CTCs "call" neutrophils (chemokines ↑ with WBCs, p = 0.001) | `ctc_cluster_wbc` samples physically contain neutrophils → CXCL8/S100A8/9 partly reflect WBC transcripts, not CTC signaling. Must de-confound. |

### ❓ Open / not yet tested
- Whether the HIF core is statistically enriched (GSEA) vs just present.
- Whether findings hold **per donor** (BR16 / LM2 / Patient) or are driven by one model.
- Whether `number of ctcs` (cluster size) correlates continuously with hypoxia score.
- Cluster identities beyond Cluster 2 (= proliferating: PRR11/CDKN3/KNL1/CDKL1).

---

## 2. The story the data currently tells

Hypoxia in CTCs is associated with the **aggressive "cluster + immune-escort" phenotype**
(clustering and neutrophil association), **not** with circadian timing. The original
"universal stress signature" and "night-time hypoxia" ideas are not supported; the
clustering/bodyguard axis is. This aligns with known CTC-cluster + neutrophil metastasis
biology and is the strongest, most positive result so far.

---

## 3. TODO (priority order, based on current results)

> **MINIMUM BAR FOR JOURNAL SUBMISSION = items 1, 2, 3 below.**
> Until these are done, the work is preprint/workshop-grade, not journal-grade.
> Items 1–3 are doable in days with existing code; items 4–5 (more data) are the hard part.


1. **[High] GSEA-validate the HIF survival core.** fgsea/Hallmark on the hypoxia-vs-rest
   overlap; report a tight enriched core, not all 847 genes. → turns "promising" into "proven."
2. **[High] De-confound the bodyguard result.** Re-test neutrophil recruitment using
   tumor-intrinsic chemokine genes only, or regress out a WBC-content score, so the signal
   is CTC-driven not WBC-contamination.
3. **[Med] Per-donor robustness.** Repeat H1b (cluster) and H2 (WBC) split by BR16/LM2/
   Patient + report effect sizes; rules out single-model artifact.
4. **[Med] Cluster-size dose-response.** Correlate hypoxia score vs `number of ctcs`
   (continuous) — strengthens the "hypoxia → clustering" claim.
5. **[Med] Annotate all 3 clusters** with markers (EMT/stemness/proliferation/hypoxia).
6. **[Low] Finalize the written draft** `paper/hypoxia_cluster_escort_paper.md`: decide
   venue, add real citations, fold in GSEA/de-confound results once done.
7. **[Low] Optional negative-result companion**: the cross-stressor benchmark can be a
   separate short workshop paper (already drafted in `paper/`).

---

## 4. Planned paper — structure

> **DRAFT WRITTEN:** `paper/hypoxia_cluster_escort_paper.md` (+ `.docx`). The structure
> below is implemented in that draft; remaining work is the analyses in §3/§5.

**Working title:** *Hypoxia Licenses the Cluster-and-Escort Phenotype of Circulating
Tumor Cells* (mechanistic framing — primary paper).

1. **Abstract** — hypoxia → CTC clustering + neutrophil escort; not circadian-timed.
2. **Introduction** — CTCs face bloodstream stress; clusters & neutrophil CTCs are most
   metastatic; question: what licenses this phenotype?
3. **Data & Methods** — GSE126669 + GSE180097; Seurat v5 + Harmony pipeline; hypoxia
   signature derivation; module scoring; statistical tests.
4. **Results**
   - 4.1 Hypoxia signature derived from GSE126669; validated within circadian data.
   - 4.2 Hypoxia is elevated in clustered CTCs (H1b). *[done; add GSEA + dose-response]*
   - 4.3 Hypoxia is elevated in neutrophil-associated CTCs (H2). *[done; add de-confound]*
   - 4.4 Shared HIF-target survival core between hypoxic and rest-phase CTCs. *[needs GSEA]*
   - 4.5 Negative control: no circadian-timing effect; stress is not universal across
     stressors (cite the cross-stressor benchmark).
5. **Discussion** — mechanistic model; clinical implication (target hypoxic CTC clusters /
   neutrophil interaction); contrast with the failed universal-stress idea.
6. **Limitations** — small hypoxia cohort (n=30); WBC composition confound addressed;
   two breast datasets only; correlational not causal.
7. **Conclusion.**

**Companion (optional, NeurIPS workshop):** *Stressor-Specific, Not Universal: a
cross-stressor transfer benchmark for CTCs* — already drafted in
`paper/cross_stressor_ctc_paper.docx` (benchmark + disentanglement GNN + negative finding).

---

## 5. Additional testing still to be done (beyond TODO)

- **GSEA / pathway enrichment** on the shared core (confirm HIF dominance).
- **WBC-content deconvolution** or regression to clean H2b.
- **Per-donor and per-cluster** stratified stats.
- **Continuous dose-response** (cluster size, WBC count) vs hypoxia score.
- **Cell-cell interaction analysis** (e.g., ligand-receptor, CellChat) between CTCs and
  attached WBCs to support "recruitment" claim mechanistically.
- **External validation** ideally: a third CTC dataset or pseudo-bulk check of the HIF core.
- **Causality caveat** stated explicitly (this is association, not perturbation).
