# ============================================
# LIBRARIES
# ============================================
library(Seurat)
library(ggplot2)
library(patchwork)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(SingleCellExperiment)

# ============================================
# HELPER FUNCTIONS
# ============================================
convert_ensembl_to_symbol <- function(expr_matrix) {
  genes   <- rownames(expr_matrix)
  symbols <- mapIds(
    org.Hs.eg.db,
    keys      = genes,
    column    = "SYMBOL",
    keytype   = "ENSEMBL",
    multiVals = "first"
  )
  symbols     <- symbols[!is.na(symbols) & symbols != ""]
  common_ids  <- intersect(genes, names(symbols))
  if (length(common_ids) == 0) stop("No Ensembl IDs could be mapped.")
  expr_sub             <- as.data.frame(expr_matrix[common_ids, , drop = FALSE])
  expr_sub$hgnc_symbol <- symbols[rownames(expr_sub)]
  expr_agg             <- aggregate(. ~ hgnc_symbol, data = expr_sub, FUN = mean)
  rownames(expr_agg)   <- expr_agg$hgnc_symbol
  expr_agg$hgnc_symbol <- NULL
  return(as.matrix(expr_agg))
}

counts_to_tpm <- function(counts_matrix, gene_lengths_kb) {
  common <- intersect(rownames(counts_matrix), names(gene_lengths_kb))
  cat("  Genes matched:", length(common), "of", nrow(counts_matrix), "\n")
  rpk <- counts_matrix[common, ] / gene_lengths_kb[common]
  tpm <- t(t(rpk) / colSums(rpk)) * 1e6
  return(tpm)
}

# ============================================
# LOAD RAW EXPRESSION - GSE126669
# ============================================
files_126669 <- list.files(
  "main/expression_data/GSE126669",
  pattern    = ".counts.txt.gz",
  full.names = TRUE
)

expr1_list <- lapply(files_126669, function(f) {
  x <- read.table(f, header = TRUE, sep = "\t")
  colnames(x) <- c("gene", basename(f))
  return(x)
})

expr1_raw                   <- Reduce(function(x, y) merge(x, y, by = "gene", all = TRUE), expr1_list)
rownames(expr1_raw)         <- expr1_raw$gene
expr1_raw$gene              <- NULL
expr1_raw[is.na(expr1_raw)] <- 0
expr1                       <- as.matrix(expr1_raw)
mode(expr1)                 <- "numeric"
colnames(expr1)             <- gsub(".counts.txt.gz", "", colnames(expr1))

rm(expr1_list, expr1_raw, files_126669)
cat("expr1 dimensions:", dim(expr1), "\n")

# ============================================
# LOAD RAW EXPRESSION - GSE180097
# ============================================
sce_br16    <- readRDS("ctc-mechanism/expression_data/GSE180097 data/GSE180097_br16.rds")
sce_lm2     <- readRDS("ctc-mechanism/expression_data/GSE180097 data/GSE180097_lm2.rds")
sce_patient <- readRDS("ctc-mechanism/expression_data/GSE180097 data/GSE180097_patient.rds")

expr3_br16    <- as.matrix(assay(sce_br16,    "counts"))
expr3_lm2     <- as.matrix(assay(sce_lm2,     "counts"))
expr3_patient <- as.matrix(assay(sce_patient, "counts"))

cat("expr1:",         dim(expr1),         "\n")
cat("expr3_br16:",    dim(expr3_br16),    "\n")
cat("expr3_lm2:",     dim(expr3_lm2),     "\n")
cat("expr3_patient:", dim(expr3_patient), "\n")

# ============================================
# ENSEMBL -> GENE SYMBOL CONVERSION
# ============================================
expr1_final         <- convert_ensembl_to_symbol(expr1)
expr3_br16_final    <- convert_ensembl_to_symbol(expr3_br16)
expr3_lm2_final     <- convert_ensembl_to_symbol(expr3_lm2)
expr3_patient_final <- convert_ensembl_to_symbol(expr3_patient)

# ============================================
# LOAD METADATA
# ============================================
meta1     <- read.csv("main/metadata/metadata_GSE126669.csv")
circ_meta <- read.csv("main/metadata/metadata_GSE180097.csv")

# Filter circ_meta to only active/resting (drops numeric timepoints like "1000","1800" etc.)
circ_meta_filtered <- circ_meta[
  circ_meta$timepoint %in% c("active", "resting"),
]
cat("circ_meta rows after filtering:", nrow(circ_meta_filtered), "\n")

# ============================================
# BUILD LABEL TABLES
# ============================================
labels_gse1 <- data.frame(
  expr_sample      = colnames(expr1_final),
  hypoxia          = meta1$hypoxia,
  model            = meta1$model,
  stringsAsFactors = FALSE
)

expr3_combined <- cbind(
  expr3_br16_final,
  expr3_lm2_final,
  expr3_patient_final
)

labels_gse3 <- data.frame(
  expr_sample      = colnames(expr3_combined),
  timepoint        = circ_meta_filtered$timepoint,
  model            = circ_meta_filtered$donor,
  stringsAsFactors = FALSE
)

cat("expr3_combined cols:", ncol(expr3_combined), "\n")
cat("labels_gse3 rows:",    nrow(labels_gse3),     "\n")

# ============================================
# LOAD POST-QC WORKSPACE
# ============================================
load("post_qc_final_workspace.RData")

# ============================================
# MANUAL LOG1P ON EXISTING SPLIT LAYERS
# NOTE: NO JoinLayers() before this step.
# PhD's instruction: FindVariableFeatures must
# see per-batch split layers, not joined data,
# otherwise HVG selection gets dominated by
# whichever batch has the strongest signal.
# ============================================
for (layer_name in Layers(seurat_merged, assay = "RNA")) {
  if (grepl("counts", layer_name)) {
    data_layer_name <- gsub("counts", "data", layer_name)
    raw_tpm          <- LayerData(seurat_merged, layer = layer_name)
    LayerData(seurat_merged, layer = data_layer_name) <- log1p(raw_tpm)
  }
}
cat("Manual log1p transformation applied to all existing split TPM layers.\n")

# ============================================
# VARIABLE FEATURES - mvp (safer for Smart-seq2 TPM)
# ============================================
seurat_merged <- FindVariableFeatures(
  seurat_merged,
  selection.method = "mvp",
  nfeatures        = 3000
)
cat("Variable features:", length(VariableFeatures(seurat_merged)), "\n")

# ============================================
# SCALE DATA (no centering - protects zero-inflated TPM)
# ============================================
seurat_merged <- ScaleData(
  seurat_merged,
  do.center = FALSE,
  do.scale  = TRUE
)
cat("Scaling completed safely for TPM.\n")

# ============================================
# PCA
# ============================================
seurat_merged <- RunPCA(seurat_merged, npcs = 50)
ElbowPlot(seurat_merged, ndims = 50)

save.image("post_pca_v2.RData")
cat("Saved post_pca_v2.RData\n")

# ============================================
# HARMONY INTEGRATION (only - no CCA)
# ============================================
cat("Running Harmony integration...\n")
seurat_harmony <- IntegrateLayers(
  object         = seurat_merged,
  method         = HarmonyIntegration,
  orig.reduction = "pca",
  new.reduction  = "harmony",
  verbose        = TRUE
)

seurat_harmony <- JoinLayers(seurat_harmony)

dims_to_use <- 1:30

cat("Running UMAP...\n")
seurat_harmony <- RunUMAP(
  seurat_harmony,
  reduction      = "harmony",
  dims           = dims_to_use,
  reduction.name = "umap.harmony"
)

cat("Finding neighbors and clusters...\n")
seurat_harmony <- FindNeighbors(
  seurat_harmony,
  reduction = "harmony",
  dims      = dims_to_use
)
seurat_harmony <- FindClusters(seurat_harmony, resolution = 0.5)
cat("Clustering complete!\n")

save.image("post_integration_v2.RData")
cat("Saved post_integration_v2.RData\n")

# ============================================
# ASSIGN LABELS
# stress_condition (binary):
#   GSE180097 -> resting = Stressed, active  = Baseline
#   GSE126669 -> hypoxia Positive   = Stressed, Negative = Baseline
# ============================================
cell_names_stripped <- gsub(
  "^(GSE126669_|BR16_|LM2_|PATIENT_)", "",
  colnames(seurat_harmony)
)

seurat_harmony$hypoxia <- labels_gse1$hypoxia[
  match(cell_names_stripped, labels_gse1$expr_sample)
]

seurat_harmony$timepoint <- labels_gse3$timepoint[
  match(cell_names_stripped, labels_gse3$expr_sample)
]

seurat_harmony$stress_condition <- ifelse(
  !is.na(seurat_harmony$timepoint),
  ifelse(seurat_harmony$timepoint == "resting", "Stressed", "Baseline"),
  ifelse(seurat_harmony$hypoxia == "Positive", "Stressed", "Baseline")
)

# stress_detail: keeps the two stressors visually distinct on one plot
seurat_harmony$stress_detail <- ifelse(
  !is.na(seurat_harmony$timepoint),
  ifelse(seurat_harmony$timepoint == "resting", "Resting (Stressed)", "Active (Baseline)"),
  ifelse(seurat_harmony$hypoxia == "Positive", "Hypoxia (Stressed)", "Normoxia (Baseline)")
)

# ============================================
# SANITY CHECKS
# ============================================
cat("\nHypoxia distribution:\n");          print(table(seurat_harmony$hypoxia,          useNA = "always"))
cat("\nTimepoint distribution:\n");        print(table(seurat_harmony$timepoint,        useNA = "always"))
cat("\nStress condition distribution:\n"); print(table(seurat_harmony$stress_condition, useNA = "always"))
cat("\nStress detail distribution:\n");    print(table(seurat_harmony$stress_detail,    useNA = "always"))

# ============================================
# PLOTS - exactly as specified by PhD
# ============================================
p1 <- DimPlot(seurat_harmony,
              reduction = "umap.harmony",
              group.by  = "dataset_id") + ggtitle("UMAP by Dataset (Integrated)")

p2 <- DimPlot(seurat_harmony,
              reduction = "umap.harmony",
              group.by  = "seurat_clusters",
              label     = TRUE) + ggtitle("UMAP by Seurat Clusters")

p1 + p2

# ============================================
# ADDITIONAL PLOTS - stress condition exploration
# ============================================

# Combined binary stress condition
p3 <- DimPlot(seurat_harmony, reduction = "umap.harmony",
              group.by = "stress_condition") + ggtitle("UMAP by Stress Condition")

# Hypoxia only (GSE126669 cells - others NA)
p4 <- DimPlot(seurat_harmony, reduction = "umap.harmony",
              group.by = "hypoxia") + ggtitle("UMAP by Hypoxia (GSE126669)")

# Timepoint only (GSE180097 cells - others NA)
p5 <- DimPlot(seurat_harmony, reduction = "umap.harmony",
              group.by = "timepoint") + ggtitle("UMAP by Timepoint (GSE180097)")

(p3 | p4) / p5

# All four conditions on one combined plot:
# resting vs active AND hypoxia vs normoxia, distinct colors
p6 <- DimPlot(
  seurat_harmony,
  reduction = "umap.harmony",
  group.by  = "stress_detail",
  cols      = c(
    "Resting (Stressed)"  = "#D55E00",
    "Active (Baseline)"   = "#56B4E9",
    "Hypoxia (Stressed)"  = "#CC79A7",
    "Normoxia (Baseline)" = "#009E73"
  )
) + ggtitle("UMAP - Resting vs Active & Hypoxia vs Normoxia")

p6

# ============================================
# SAVE FINAL WORKSPACE
# ============================================
save.image("post_integration_stress_v2.RData")
cat("Saved post_integration_stress_v2.RData\n")
