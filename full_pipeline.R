# ============================================
# LIBRARIES
# ============================================
library(Seurat)
library(ggplot2)
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
  expr_sub            <- as.data.frame(expr_matrix[common_ids, , drop = FALSE])
  expr_sub$hgnc_symbol <- symbols[rownames(expr_sub)]
  expr_agg            <- aggregate(. ~ hgnc_symbol, data = expr_sub, FUN = mean)
  rownames(expr_agg)  <- expr_agg$hgnc_symbol
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
  pattern = ".counts.txt.gz",
  full.names = TRUE
)

expr1_list <- lapply(files_126669, function(f) {
  x <- read.table(f, header = TRUE, sep = "\t")
  colnames(x) <- c("gene", basename(f))
  return(x)
})

expr1_raw             <- Reduce(function(x, y) merge(x, y, by = "gene", all = TRUE), expr1_list)
rownames(expr1_raw)   <- expr1_raw$gene
expr1_raw$gene        <- NULL
expr1_raw[is.na(expr1_raw)] <- 0
expr1                 <- as.matrix(expr1_raw)
mode(expr1)           <- "numeric"
colnames(expr1)       <- gsub(".counts.txt.gz", "", colnames(expr1))

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

cat("expr1:",        dim(expr1),        "\n")
cat("expr3_br16:",   dim(expr3_br16),   "\n")
cat("expr3_lm2:",    dim(expr3_lm2),    "\n")
cat("expr3_patient:", dim(expr3_patient), "\n")

# ============================================
# ENSEMBL -> GENE SYMBOL CONVERSION
# ============================================
expr1_final        <- convert_ensembl_to_symbol(expr1)
expr3_br16_final   <- convert_ensembl_to_symbol(expr3_br16)
expr3_lm2_final    <- convert_ensembl_to_symbol(expr3_lm2)
expr3_patient_final <- convert_ensembl_to_symbol(expr3_patient)

# ============================================
# LOAD METADATA
# ============================================
meta1     <- read.csv("main/metadata/metadata_GSE126669.csv")
circ_meta <- read.csv("main/metadata/metadata_GSE180097.csv")

# Filter circ_meta to only active/resting
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
cat("labels_gse3 rows:",   nrow(labels_gse3),    "\n")

# ============================================
# LOAD POST-QC WORKSPACE AND BUILD SEURAT
# ============================================
load("post_qc_final_workspace.RData")

# ============================================
# PRE-PROCESSING
# ============================================
seurat_merged <- JoinLayers(seurat_merged)

seurat_merged <- SetAssayData(
  seurat_merged,
  layer    = "data",
  new.data = log1p(GetAssayData(seurat_merged, layer = "counts"))
)
cat("Log1p done\n")

seurat_merged[["RNA"]] <- split(
  seurat_merged[["RNA"]],
  f = seurat_merged$dataset_id
)
cat("Layers:", paste(Layers(seurat_merged), collapse = ", "), "\n")

seurat_merged <- FindVariableFeatures(seurat_merged, nfeatures = 3000)
cat("Variable features:", length(VariableFeatures(seurat_merged)), "\n")

seurat_merged <- ScaleData(seurat_merged)
cat("Scaling done\n")

seurat_merged <- RunPCA(seurat_merged, npcs = 50)
ElbowPlot(seurat_merged, ndims = 50)

save.image("post_pca.RData")
cat("Saved post_pca.RData\n")

# ============================================
# HARMONY INTEGRATION
# ============================================
seurat_harmony <- IntegrateLayers(
  object         = seurat_merged,
  method         = HarmonyIntegration,
  orig.reduction = "pca",
  new.reduction  = "harmony",
  verbose        = TRUE
)

seurat_harmony <- JoinLayers(seurat_harmony)
seurat_harmony <- FindNeighbors(seurat_harmony, reduction = "harmony", dims = 1:20)
seurat_harmony <- FindClusters(seurat_harmony, resolution = 0.5)
seurat_harmony <- RunUMAP(seurat_harmony, reduction = "harmony", dims = 1:20)

# ============================================
# ASSIGN LABELS - HARMONY
# ============================================
cell_names_stripped_harmony <- gsub(
  "^(GSE126669_|BR16_|LM2_|PATIENT_)", "",
  colnames(seurat_harmony)
)

seurat_harmony$hypoxia <- labels_gse1$hypoxia[
  match(cell_names_stripped_harmony, labels_gse1$expr_sample)
]

seurat_harmony$timepoint <- labels_gse3$timepoint[
  match(cell_names_stripped_harmony, labels_gse3$expr_sample)
]

# stress_condition: timepoint for GSE180097, hypoxia for GSE126669
seurat_harmony$stress_condition <- seurat_harmony$timepoint
seurat_harmony$stress_condition[
  is.na(seurat_harmony$stress_condition)
] <- ifelse(
  seurat_harmony$hypoxia[
    is.na(seurat_harmony$stress_condition)
  ] == "Positive",
  "Hypoxia",
  "Normoxia"
)

cat("\nHarmony - Hypoxia distribution:\n")
print(table(seurat_harmony$hypoxia,         useNA = "always"))
cat("\nHarmony - Timepoint distribution:\n")
print(table(seurat_harmony$timepoint,        useNA = "always"))
cat("\nHarmony - Stress condition distribution:\n")
print(table(seurat_harmony$stress_condition, useNA = "always"))

# ============================================
# DIMPLOTS - HARMONY
# ============================================
p1 <- DimPlot(seurat_harmony, reduction = "umap",
              group.by = "dataset_id")        + ggtitle("Harmony - by Dataset")
p2 <- DimPlot(seurat_harmony, reduction = "umap",
              group.by = "seurat_clusters",
              label    = TRUE)                + ggtitle("Harmony - by Cluster")
p3 <- DimPlot(seurat_harmony, reduction = "umap",
              group.by = "stress_condition")  + ggtitle("Harmony - Stress Condition")
(p1 | p2) / p3

save.image("after_stress_after_harmony_final.RData")
cat("Saved harmony workspace\n")

# ============================================
# CCA INTEGRATION
# (load from post_pca to keep separate)
# ============================================
load("post_pca.RData")

seurat_merged <- JoinLayers(seurat_merged)

seurat_merged <- SetAssayData(
  seurat_merged,
  layer    = "data",
  new.data = log1p(GetAssayData(seurat_merged, layer = "counts"))
)

seurat_merged[["RNA"]] <- split(
  seurat_merged[["RNA"]],
  f = seurat_merged$dataset_id
)

cat("Layers after split:\n")
print(Layers(seurat_merged))

seurat_cca <- IntegrateLayers(
  object         = seurat_merged,
  method         = CCAIntegration,
  orig.reduction = "pca",
  new.reduction  = "integrated.cca",
  dims           = 1:20,
  k.anchor       = 3,
  k.filter       = 10,
  k.score        = 10,
  k.weight       = 20,
  verbose        = TRUE
)

seurat_cca <- JoinLayers(seurat_cca)
cat("CCA integration done\n")

seurat_cca <- FindNeighbors(seurat_cca, reduction = "integrated.cca", dims = 1:20)
seurat_cca <- FindClusters(seurat_cca,  resolution = 0.5)
seurat_cca <- RunUMAP(
  seurat_cca,
  reduction      = "integrated.cca",
  dims           = 1:20,
  reduction.name = "umap.cca"
)
cat("UMAP done\n")

# ============================================
# ASSIGN LABELS - CCA
# ============================================
cell_names_stripped_cca <- gsub(
  "^(GSE126669_|BR16_|LM2_|PATIENT_)", "",
  colnames(seurat_cca)
)

seurat_cca$hypoxia <- labels_gse1$hypoxia[
  match(cell_names_stripped_cca, labels_gse1$expr_sample)
]

seurat_cca$timepoint <- labels_gse3$timepoint[
  match(cell_names_stripped_cca, labels_gse3$expr_sample)
]

# stress_condition: timepoint for GSE180097, hypoxia for GSE126669
seurat_cca$stress_condition <- seurat_cca$timepoint
seurat_cca$stress_condition[
  is.na(seurat_cca$stress_condition)
] <- ifelse(
  seurat_cca$hypoxia[
    is.na(seurat_cca$stress_condition)
  ] == "Positive",
  "Hypoxia",
  "Normoxia"
)

cat("\nCCA - Hypoxia distribution:\n")
print(table(seurat_cca$hypoxia,         useNA = "always"))
cat("\nCCA - Timepoint distribution:\n")
print(table(seurat_cca$timepoint,        useNA = "always"))
cat("\nCCA - Stress condition distribution:\n")
print(table(seurat_cca$stress_condition, useNA = "always"))

# ============================================
# DIMPLOTS - CCA
# ============================================
p1 <- DimPlot(seurat_cca, reduction = "umap.cca",
              group.by = "dataset_id")        + ggtitle("CCA - by Dataset")
p2 <- DimPlot(seurat_cca, reduction = "umap.cca",
              group.by = "seurat_clusters",
              label    = TRUE)                + ggtitle("CCA - by Cluster")
p3 <- DimPlot(seurat_cca, reduction = "umap.cca",
              group.by = "hypoxia")           + ggtitle("CCA - Hypoxia (GSE126669)")
p4 <- DimPlot(seurat_cca, reduction = "umap.cca",
              group.by = "timepoint")         + ggtitle("CCA - Timepoint (GSE180097)")
p5 <- DimPlot(seurat_cca, reduction = "umap.cca",
              group.by = "stress_condition")  + ggtitle("CCA - Stress Condition")
(p1 | p2) / (p3 | p4) / p5

save.image("post_cca_integration.RData")
cat("Saved post_cca_integration.RData\n")