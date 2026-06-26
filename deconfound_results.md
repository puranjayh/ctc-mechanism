# De-confounding + GSEA results (2026-06-24)

Run from `ctc_obj.rds` via `gsea_deconfound.R` and `perdonor_breakdown.R`.
These results **materially revise** the project's original mechanistic hypotheses.

## A. GSEA validation of the hypoxia survival core — HOLDS ✅

Ranked all 10,558 detected genes by the hypoxia (Positive vs Negative, GSE126669)
contrast and ran fgsea vs MSigDB Hallmark.

| Set | NES | p | padj |
|---|---|---|---|
| HIF_CORE_HANDOFF (curated 9-gene) | 1.68 | 0.009 | **0.044** |
| HALLMARK_HYPOXIA | 1.57 | 0.002 | **0.021** |

- The curated HIF core is the **#1 enriched Hallmark-level set**.
- Tight, defensible core = the 43-gene **HALLMARK_HYPOXIA leading edge**
  (`hif_core_leadingedge.csv`) — use this instead of the loose 847-gene list.
  Recovers VEGFA/P4HA1/NDRG1/HK2/ERO1A; HIF1A/EGLN1/BNIP3/MTHFD2 are simply
  not members of MSigDB Hallmark hypoxia.
- This signal lives **inside the hypoxia dataset**, so it is NOT affected by the
  donor confound below. It is the one solid positive result.

## B. H2b "neutrophil recruitment" — CONTAMINATION, drop the claim ❌

- Of the 11 genes in the original `Neutrophil1` recruitment score, **only S100A9
  survives QC** in this sparse Smart-seq2 object (cor(score, S100A9) = 0.993).
- Tumor-secreted recruiters needed to test the hypothesis — CSF3, CSF2,
  CXCL1/2/3/5/6/8 — are **not detected at all**.
- S100A9 is a neutrophil-INTRINSIC gene; it is high in `ctc_cluster_wbc` because
  those samples physically contain neutrophils. This is contamination, not CTC
  signaling. The hypothesis is **untestable** with this gene panel.

## C. H1b (cluster) and H2a (WBC) — DONOR ARTIFACTS (Simpson's paradox) ❌

The three sources have very different baseline hypoxia scores, and cluster/WBC
cells are unevenly distributed across them:

```
                BR16    LM2   PATIENT
baseline hyp   -0.013  -0.333  +0.392
with_WBC  n      43      9       6     (74% of WBC cells are BR16)
no_WBC    n      94     97      27
```

**H2a per donor (hypoxia, with_WBC vs no_WBC):**
```
BR16     -0.000 vs -0.019   p = 0.62
LM2      -0.335 vs -0.333   p = 0.86
PATIENT  +0.393 vs +0.392   p = 0.84
```
Pooled p = 0.008 → within every donor the effect is ZERO. With donor controlled:
`lm` wbc coef = 0.013, p = 0.638. **Not real.**

**H1b per donor (hypoxia, cluster vs single):**
```
BR16     p = 0.57
LM2      p = 0.046  (cluster LOWER — wrong direction)
PATIENT  p = 0.64
```
With donor controlled: cluster coef = -0.013, p = 0.558. **Not real.**

(`WBCcontent` adjustment alone did NOT kill H2a — p stayed 0.008 — it was the
**donor** term that did. cor(HypoxiaSig, WBCcontent) = 0.20, so no circularity.)

Figure: `deconfound_simpson.png` (pooled illusion vs within-donor truth).

## Revised scoreboard

| Claim | Old verdict | After de-confounding |
|---|---|---|
| Hypoxia survival core (GSEA) | core proposed | ✅ **validated** (padj 0.02) |
| H1 hypoxia ↑ at rest | ❌ p=0.41 | ❌ (unchanged) |
| H1b hypoxia ↑ in clusters | ✅ p=0.023 | ❌ **donor artifact** |
| H2a hypoxia ↑ with WBCs | ✅ p=0.008 | ❌ **donor artifact** |
| H2b chemokines ↑ with WBCs | ⚠️ p=0.001 | ❌ **contamination (S100A9 only)** |

## D. Metastasis timing (advisor request) — a NEW robust positive ✅

`metastasis_timing.R`. GSE180097 is the Aceto sleep-metastasis dataset; rest is
the pro-metastatic window. Question: do metastatic phenotypes peak at rest AFTER
donor control? (proliferation panel = 22 cell-cycle genes.)

**Proliferation ↑ at rest — robust, unlike the hypoxia results:**
```
              rest    active   p
pooled        0.165   -0.355   0.015
+ donor                        0.0036   <- STRENGTHENS with donor control
BR16          0.648   -0.104   0.022    (sig)
LM2          -0.080   -0.486   0.47     (same direction)
PATIENT      -0.711   -1.249   0.70     (same direction)
```
Same direction in ALL three donors; donor control makes it stronger (the opposite
of Simpson's paradox). Reproduces the published sleep -> proliferation biology.

- **Clustering ↑ at rest:** 66% vs 57% clustered; consistent across donors but
  only a trend (donor-controlled logistic p = 0.077; BR16 fisher p = 0.058).
- **Hypoxia vs phase:** flat (p = 0.41); with donor control trends LOWER at rest
  (p = 0.059) — opposite to the advisor's hypoxia<->rest bridge. Confirmed NOT involved.

Caveat: proliferation panel included 4 cluster-2 markers (PRR11/CDKN3/KNL1/CDKL1)
from the handoff alongside 18 canonical cell-cycle genes; result is driven by the
canonical set, but a cell-cycle-only re-run is a cheap robustness check.

## Bottom line / the honest, de-confounded story

The original "hypoxia -> clustering + neutrophil escort" model does **not** survive.
But a cleaner, properly-controlled story emerges:

1. **Rest-phase CTCs are more proliferative** (donor-robust, p=0.0036) and trend
   toward more clustered — i.e. the metastatic-timing signal is real and
   reproduces known biology. [NEW positive]
2. This is **hypoxia-INDEPENDENT** — hypoxia does not track rest phase, clustering,
   or WBC attachment once donor is controlled (all were Simpson's-paradox artifacts).
3. There IS a bona-fide **hypoxia survival core, GSEA-validated** (HALLMARK_HYPOXIA
   padj=0.02) — but it operates within the hypoxia axis, not the circadian/cluster axis.
4. Neutrophil recruitment is **untestable** here (signal = S100A9 contamination).
5. Cross-stressor transfer remains **negative** (settled).

So: two independent, defensible findings (rest->proliferation timing; a validated
hypoxia survival program) plus a rigorous de-confounding that separates real signal
from cell-line/contamination artifacts — a stronger and more honest paper than the
original mechanistic claim.
