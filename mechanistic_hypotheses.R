# ============================================================
# MECHANISTIC HYPOTHESES (PhD direction)
#   Model: hypoxia DRIVES rest-phase shedding + WBC recruitment
#   H1 Night-Time Hypoxia : hypoxia program is higher in
#        resting (vs active) and in clustered (vs single) CTCs
#   H2 Bodyguard          : hypoxic CTCs recruit neutrophils
#        (higher hypoxia score + chemokines when WBCs attached)
# ============================================================
library(Seurat); library(ggplot2)
load("post_integration_stress_v2.RData")
obj <- JoinLayers(seurat_harmony)

# ------------------------------------------------------------
# 1. ATTACH GSE180097 metadata (phase, sample type, clumping, WBC)
#    (positional map: expr3_combined cols <-> circ_meta_filtered rows,
#     same convention used in your pipeline; cross-tab below verifies)
# ------------------------------------------------------------
circ <- read.csv("main/metadata/metadata_GSE180097.csv")
circ_f <- circ[circ$timepoint %in% c("active", "resting"), ]
if (!exists("expr3_combined"))
  expr3_combined <- cbind(expr3_br16_final, expr3_lm2_final, expr3_patient_final)

lab <- data.frame(
  expr_sample = colnames(expr3_combined),
  phase       = circ_f$timepoint,
  sample_type = circ_f$sample.type,
  donor       = circ_f$donor,
  n_ctc       = suppressWarnings(as.numeric(circ_f$number.of.ctcs)),
  n_wbc       = suppressWarnings(as.numeric(circ_f$number.of.wbcs.attached)),
  stringsAsFactors = FALSE
)
bc <- gsub("^(GSE126669_|BR16_|LM2_|PATIENT_)", "", colnames(obj))
m  <- match(bc, lab$expr_sample)
obj$phase       <- lab$phase[m]                                   # resting / active
obj$sample_type <- lab$sample_type[m]                             # ctc_single / ctc_cluster / ctc_cluster_wbc
obj$clumping    <- ifelse(grepl("cluster", obj$sample_type), "cluster", "single")
obj$wbc_status  <- ifelse(obj$sample_type == "ctc_cluster_wbc", "with_WBC",
                   ifelse(!is.na(obj$sample_type), "no_WBC", NA))
cat("Sample types mapped:\n"); print(table(obj$sample_type, useNA = "always"))
cat("Phase x WBC:\n");        print(table(obj$phase, obj$wbc_status, useNA = "always"))

# ------------------------------------------------------------
# 2. DATA-DRIVEN HYPOXIA SURVIVAL SIGNATURE (from GSE126669)
#    -> top genes up in hypoxia-Positive single CTCs
# ------------------------------------------------------------
hyp <- subset(obj, cells = colnames(obj)[obj$hypoxia %in% c("Positive", "Negative")])
Idents(hyp) <- hyp$hypoxia
hyp_mk <- FindMarkers(hyp, ident.1 = "Positive", ident.2 = "Negative",
                      only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
hyp_sig <- head(rownames(hyp_mk[order(-hyp_mk$avg_log2FC), ]), 50)
write.csv(hyp_mk, "hypoxia_survival_markers.csv")
cat("\nHypoxia survival signature (top genes):\n"); print(head(hyp_sig, 20))

obj <- AddModuleScore(obj, features = list(intersect(hyp_sig, rownames(obj))),
                      name = "HypoxiaSig", ctrl = 50)

# ------------------------------------------------------------
# 3. H1  NIGHT-TIME HYPOXIA  (within GSE180097 only)
# ------------------------------------------------------------
circ_cells <- subset(obj, cells = colnames(obj)[!is.na(obj$phase)])
cat("\n== H1: hypoxia score, RESTING vs ACTIVE ==\n")
print(tapply(circ_cells$HypoxiaSig1, circ_cells$phase, mean))
print(wilcox.test(HypoxiaSig1 ~ phase, data = circ_cells@meta.data))
cat("\n== H1b: hypoxia score, CLUSTER vs SINGLE ==\n")
print(tapply(circ_cells$HypoxiaSig1, circ_cells$clumping, mean))
print(wilcox.test(HypoxiaSig1 ~ clumping, data = circ_cells@meta.data))

p_h1 <- VlnPlot(circ_cells, "HypoxiaSig1", group.by = "phase", pt.size = 0.3) +
        ggtitle("H1: Hypoxia signature by phase") + NoLegend()
p_veg <- VlnPlot(circ_cells, "VEGFA", group.by = "phase", pt.size = 0.3) +
        ggtitle("VEGFA by phase") + NoLegend()

# ------------------------------------------------------------
# 4. H2  BODYGUARD  (hypoxic CTCs recruit neutrophils)
# ------------------------------------------------------------
neutro <- c("CXCL8","CXCL1","CXCL2","CXCL3","CXCL5","CXCL6","CSF3","CSF2",
            "S100A8","S100A9","IL6")
obj <- AddModuleScore(obj, features = list(intersect(neutro, rownames(obj))),
                      name = "Neutrophil", ctrl = 50)
circ_cells <- subset(obj, cells = colnames(obj)[!is.na(obj$phase)])
wbc_cells  <- subset(circ_cells, cells = colnames(circ_cells)[!is.na(circ_cells$wbc_status)])

cat("\n== H2: hypoxia score, WITH_WBC vs NO_WBC ==\n")
print(tapply(wbc_cells$HypoxiaSig1, wbc_cells$wbc_status, mean))
print(wilcox.test(HypoxiaSig1 ~ wbc_status, data = wbc_cells@meta.data))
cat("\n== H2b: neutrophil-chemokine score, WITH_WBC vs NO_WBC ==\n")
print(tapply(wbc_cells$Neutrophil1, wbc_cells$wbc_status, mean))
print(wilcox.test(Neutrophil1 ~ wbc_status, data = wbc_cells@meta.data))

p_h2  <- VlnPlot(wbc_cells, "HypoxiaSig1", group.by = "wbc_status", pt.size = 0.3) +
         ggtitle("H2: Hypoxia score by WBC attachment") + NoLegend()
p_h2b <- VlnPlot(wbc_cells, "Neutrophil1", group.by = "wbc_status", pt.size = 0.3) +
         ggtitle("H2: Neutrophil chemokines by WBC attachment") + NoLegend()

# ------------------------------------------------------------
# 5. CLUSTER MARKERS  (is each cluster hypoxia- or phase-defined?)
# ------------------------------------------------------------
Idents(obj) <- obj$seurat_clusters
cluster_markers <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(cluster_markers, "cluster_markers.csv")
cat("\nTop 5 markers per cluster:\n")
print(do.call(rbind, by(cluster_markers, cluster_markers$cluster,
      function(d) head(d[order(-d$avg_log2FC), c("cluster","gene","avg_log2FC")], 5))))

# ------------------------------------------------------------
# 6. SHARED SURVIVAL GENES: hypoxia-Positive vs rest-phase CTCs
# ------------------------------------------------------------
rest <- subset(circ_cells, cells = colnames(circ_cells)[circ_cells$phase %in% c("resting","active")])
Idents(rest) <- rest$phase
rest_mk <- FindMarkers(rest, ident.1 = "resting", ident.2 = "active",
                       only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
shared <- intersect(rownames(hyp_mk), rownames(rest_mk))
cat("\n== Shared survival genes (hypoxia-Positive AND resting-up) ==\n")
print(shared)
write.csv(data.frame(gene = shared), "shared_survival_genes.csv", row.names = FALSE)

# ------------------------------------------------------------
# DISPLAY + SAVE
# ------------------------------------------------------------
print((p_h1 | p_veg) / (p_h2 | p_h2b))
save.image("post_mechanistic_hypotheses.RData")
cat("\nSaved -> post_mechanistic_hypotheses.RData\n")
cat("CSVs: hypoxia_survival_markers.csv, cluster_markers.csv, shared_survival_genes.csv\n")
# ------------------------------------------------------------
# READ:
#  H1 holds if hypoxia score is higher in RESTING (p<0.05) and/or CLUSTER.
#  H2 holds if hypoxia score AND neutrophil chemokines are higher in WITH_WBC.
#  'shared_survival_genes.csv' = genes used by BOTH hypoxic and rest-phase CTCs.
# ------------------------------------------------------------
