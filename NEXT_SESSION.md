# Next Session — What To Do (feed this to Claude on the new laptop)

> Paste this whole file into a new Claude session to resume work. Pair with
> `PROJECT_HANDOFF.md` (full context) and `current_state.md` (results + TODO).
> Written 2026-06-19.

---

## Short answer: yes, do items 1–3 now

These three analyses are **safe to run regardless of what my PhD decides** about
framing/venue — any version of the paper needs them, and reviewers will demand them.
They turn "promising leads" into "validated findings" and move the work from
preprint/workshop-grade to journal-grade. So start them now; don't sit idle waiting.

**But first:** confirm `ctc_obj.rds` loads cleanly on this laptop before any heavy work.
```r
obj <- readRDS("ctc_obj.rds")   # verify the data transferred
```
Also make sure the `expression_data/` folder + metadata transferred (gitignored, came via cloud/USB).

---

## What to wait for my PhD on (vs do now)

- **Wait for her on:** the direction/framing — which story to lead with, which venue,
  whether to push the circadian-timing angle further. Those are her calls.
- **Do NOT wait on:** items 1–3 below. They strengthen ANY framing. Start with item 2
  (de-confounding) — it's the biggest weakness and could sink the bodyguard claim.

---

## The three journal-minimum tasks (code-wise, existing data only)

### 1. GSEA-validate the HIF survival core  [HIGH]
- Input: `shared_survival_genes.csv` (and/or `hypoxia_survival_markers.csv`).
- Run `fgsea` against MSigDB Hallmark gene sets.
- Goal: show `HALLMARK_HYPOXIA` (and related) is **statistically enriched**, not just
  "these genes appear." Report the tight enriched core, not all 847 genes.

### 2. De-confound the bodyguard result  [HIGHEST PRIORITY]
- Problem: `ctc_cluster_wbc` samples physically contain neutrophils, so the
  neutrophil-chemokine signal (p=0.001) partly reflects WBC transcripts, not CTC signaling.
- Fix (either/both): (a) re-test recruitment using ONLY tumor-intrinsic genes, excluding
  neutrophil/myeloid genes the WBCs contribute; and/or (b) compute a "WBC-content score"
  per cell and regress it out / use as covariate, then re-test the hypoxia–WBC association.
- Goal: show the hypoxia–WBC link (p=0.008) survives, and decide whether any CTC-driven
  recruitment signal remains after removing WBC contamination.
- Modifies `mechanistic_hypotheses.R`.

### 3. Per-donor robustness  [MED]
- Repeat the clumping test (hypoxia: cluster vs single, was p=0.023) and the WBC test
  (was p=0.008) SEPARATELY for BR16, LM2, and Patient.
- Report effect sizes per donor. Goal: prove it's not one cell line driving everything.

---

## Nice-to-have (after 1–3, strengthens further)
4. Cluster-size dose-response: correlate hypoxia score vs `number of ctcs` (continuous).
5. Annotate all 3 Louvain clusters (EMT/stemness/proliferation/hypoxia markers).
6. Ligand–receptor analysis (CellChat) between CTCs and attached WBCs — makes the
   "recruitment" claim mechanistic.

## Hard / needs more data (not code-only)
- External validation in a third CTC dataset.
- The n=30 hypoxia-cohort limitation (signature rests on 30 cells).

---

## Settled facts — do NOT redo
- Cross-stressor transfer is a SETTLED NEGATIVE (4 methods ≈ chance). Don't re-run or
  tune the GNN — that would be p-hacking.
- Circadian timing does NOT explain hypoxia (p=0.41). Not supported.

## Suggested opening instruction for the new session
> "Read PROJECT_HANDOFF.md, current_state.md, and NEXT_SESSION.md. I want to do TODO
>  items 1–3 (journal-minimum). Start with item 2 (de-confound the bodyguard result).
>  Data is in ctc_obj.rds. Write the script, I'll run it and paste the output."
