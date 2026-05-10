# ctc-mechanism
# CTC Multi-Stress Immune Evasion Framework
Integrated transcriptomic analysis of circulating tumor cells 
across hypoxic, immune, and circadian stress conditions to identify 
pan-stress survival hub genes and predict immune evasion mechanisms 
using graph neural networks.

## What This Does
- Integrates three scRNA-seq CTC datasets from the Aceto Lab, ETH Zurich
- Performs batch correction across datasets using Harmony
- Identifies genes consistently upregulated across all three stress conditions
- Builds co-expression network via WGCNA to identify hub genes
- Trains Random Forest and GNN classifiers to predict stress condition from gene expression
- Maps predicted stress condition to immune evasion mechanism using literature-derived signatures
- Validates hub genes as prognostic markers in TCGA breast cancer cohort

## Datasets
| Dataset | Journal | Stress Condition | Samples |
|---------|---------|-----------------|---------|
| GSE126669 | Cell Reports 2020 | Hypoxic vs normoxic CTC clusters | 31 |
| GSE109761 | Nature 2019 | Single CTCs vs clusters under immune pressure | 357 |
| GSE180097 | Nature 2022 | Resting vs active phase CTCs (sleep stress) | 306 |

## Key Results
[To be updated upon completion]

## Biological Background
CTCs face three distinct survival stresses during hematogenous dissemination.
Existing studies examine each stress independently — no framework has integrated
transcriptomic data across conditions to identify conserved survival mechanisms.
Full literature review covers CCL5/Treg evasion axis (Sun et al. 2021),
surfaceome-based immune escape (MHC-I downregulation, aberrant glycosylation),
and hypoxia-driven CTC cluster formation (Donato et al. 2020).

## Pipeline
