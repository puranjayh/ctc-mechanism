# ============================================================
# EXPORT for the GNN (run in R/RStudio)
# Dumps 3 small CSVs the Python script reads:
#   node_features.csv  - Harmony embedding per cell (batch-corrected)
#   edges.csv          - cell-cell graph (Seurat SNN) as edge list
#   node_meta.csv      - labels: dataset, stress, hypoxia, timepoint
# ============================================================
library(Seurat)
load("post_integration_stress_v2.RData")
obj <- seurat_harmony

# ---- node features: Harmony embedding (low-dim, batch-corrected) ----
emb <- Embeddings(obj, reduction = "harmony")          # cells x dims
write.csv(data.frame(barcode = rownames(emb), emb),
          "node_features.csv", row.names = FALSE)

# ---- graph: pick the SNN graph (fallback to first available) ----
gname <- grep("snn", names(obj@graphs), value = TRUE, ignore.case = TRUE)[1]
if (is.na(gname)) gname <- names(obj@graphs)[1]
cat("Using graph:", gname, "\n")
G  <- obj@graphs[[gname]]
tt <- as(G, "TsparseMatrix")                            # i, j, x triplets
edges <- data.frame(source = rownames(G)[tt@i + 1],
                    target = colnames(G)[tt@j + 1],
                    weight = tt@x)
edges <- edges[edges$weight > 0 & edges$source != edges$target, ]
write.csv(edges, "edges.csv", row.names = FALSE)

# ---- metadata / labels ----
meta <- data.frame(
  barcode    = colnames(obj),
  dataset_id = obj$dataset_id,
  stress     = obj$stress_condition,           # Stressed / Baseline
  hypoxia    = obj$hypoxia,                     # Positive/Negative or NA
  timepoint  = obj$timepoint                    # resting/active or NA
)
write.csv(meta, "node_meta.csv", row.names = FALSE)

cat(sprintf("Exported: %d cells, %d edges, %d feature dims\n",
            nrow(emb), nrow(edges), ncol(emb)))
