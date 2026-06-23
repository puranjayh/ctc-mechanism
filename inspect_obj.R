# Quick structural inspection of ctc_obj.rds
suppressMessages(library(Seurat))
obj <- readRDS("ctc_obj.rds")

cat("=== Object summary ===\n"); print(obj)
cat("\nN genes:", nrow(obj), "  N cells:", ncol(obj), "\n")
cat("DefaultAssay:", DefaultAssay(obj), "\n")
cat("\n=== Assays ===\n"); print(Assays(obj))
cat("\n=== Layers (default assay) ===\n"); print(Layers(obj))
cat("\n=== Reductions ===\n"); print(Reductions(obj))

cat("\n=== meta.data columns ===\n"); print(colnames(obj@meta.data))

cat("\n=== key metadata tables ===\n")
for (col in c("dataset_id","orig.ident","hypoxia","phase","sample_type",
              "clumping","wbc_status","seurat_clusters","donor","stress_condition")) {
  if (col %in% colnames(obj@meta.data)) {
    cat("\n--", col, "--\n"); print(table(obj@meta.data[[col]], useNA = "always"))
  } else cat("\n--", col, ": ABSENT --\n")
}

cat("\n=== module-score columns present? ===\n")
print(grep("Sig|Neutro|Score|Module", colnames(obj@meta.data), value = TRUE))

cat("\n=== HIF core genes present in rownames? ===\n")
hif <- c("VEGFA","P4HA1","NDRG1","HK2","HIF1A","EGLN1","BNIP3","ERO1A","MTHFD2")
print(setNames(hif %in% rownames(obj), hif))

cat("\n=== data layer sanity (is it log-normalized?) ===\n")
d <- tryCatch(as.matrix(LayerData(obj, "data")[1:min(50, nrow(obj)), 1:min(5, ncol(obj))]),
              error = function(e) NULL)
if (!is.null(d)) cat("range of a data slice:", range(d), "\n") else cat("no 'data' layer accessible\n")
