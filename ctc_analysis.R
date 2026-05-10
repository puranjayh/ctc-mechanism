# ============================================
# CTC Multi-Stress Transcriptomic Analysis
# Author: Puranjay Haldankar
# Supervised by: Dr. Sarita Poonia
# ============================================

# ---- 1. Load Libraries ----
library(Seurat)
library(SingleCellExperiment)
library(harmony)
library(dplyr)
library(ggplot2)

cat("Libraries loaded\n")

# ---- 2. Set Working Directory ----
setwd("/cloud/project/ctc-mechanism")

# ---- 3. Load Expression Data ----

# GSE126669 - Hypoxic vs Normoxic CTC clusters
expr1 <- read.table(
  "expression_data/GSE126669_tpm_length_scaled_matrix.txt.gz",
  header=TRUE, row.names=1, sep="\t"
)
cat("GSE126669:", dim(expr1), "\n")

# GSE109761 - Single vs Cluster CTCs
expr2 <- read.table(
  "expression_data/GSE109761_processed_normalized_matrix_hs.txt.gz",
  header=TRUE, row.names=1, sep="\t"
)
cat("GSE109761:", dim(expr2), "\n")

# GSE180097 - Resting vs Active CTCs
sce_patient <- readRDS("expression_data/GSE180097 data/GSE180097_patient.rds")
sce_br16 <- readRDS("expression_data/GSE180097 data/GSE180097_br16.rds")
sce_lm2 <- readRDS("expression_data/GSE180097 data/GSE180097_lm2.rds")

cat("GSE180097 patient:", dim(sce_patient), "\n")
cat("GSE180097 br16:", dim(sce_br16), "\n")
cat("GSE180097 lm2:", dim(sce_lm2), "\n")

# ---- 4. Load Metadata ----
meta1 <- read.csv("metadata/metadata_GSE126669.csv")
meta2 <- read.csv("metadata/metadata_GSE109761.csv")
meta3 <- read.csv("metadata/metadata_GSE180097.csv")

cat("Metadata loaded\n")
cat("GSE126669 metadata columns:", colnames(meta1), "\n")
cat("GSE109761 metadata columns:", colnames(meta2), "\n")
cat("GSE180097 metadata columns:", colnames(meta3), "\n")
