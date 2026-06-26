# ============================================================
#  GSEA VALIDATION + H2b DE-CONFOUNDING
#  Loads ctc_obj.rds (NOT the broken .RData).
#
#  PART A. GSEA-validate the hypoxia survival core
#          - rank ALL genes by the hypoxia (Pos vs Neg) contrast
#          - fgsea vs MSigDB Hallmark  -> expect HALLMARK_HYPOXIA up
#          - report the leading-edge core (replaces the loose 847-gene list)
#
#  PART B. De-confound the neutrophil-chemokine signal (H2b)
#          The 'ctc_cluster_wbc' samples physically CONTAIN neutrophils,
#          so the raw chemokine signal is partly WBC transcript, not CTC
#          "calling out". We separate:
#            - WBCcontent  : leukocyte IDENTITY markers (the contaminant)
#            - NeutroAuto  : chemokines neutrophils express themselves
#            - NeutroTumor : chemokines a TUMOR cell would secrete to recruit
#          then ask whether (i) the with_WBC effect and (ii) the HYPOXIA
#          effect on tumor-intrinsic recruitment survive controlling for
#          WBC content.
# ============================================================

suppressMessages({
  library(Seurat); library(ggplot2); library(fgsea); library(msigdbr); library(dplyr)
})
set.seed(1)

obj <- readRDS("ctc_obj.rds")
DefaultAssay(obj) <- "RNA"
stopifnot(all(c("hypoxia","wbc_status","HypoxiaSig1") %in% colnames(obj@meta.data)))

# HIF survival core proposed in the handoff (the thing we want to validate)
hif_core <- c("VEGFA","P4HA1","NDRG1","HK2","HIF1A","EGLN1","BNIP3","ERO1A","MTHFD2")

# ============================================================
#  PART A — GSEA
# ============================================================
cat("\n##################  PART A : GSEA  ##################\n")

# --- A1. sanity: does the HypoxiaSig module score separate Pos vs Neg? ---
hyp_cells <- obj@meta.data[obj$hypoxia %in% c("Positive","Negative"), ]
cat("\n[A1] HypoxiaSig1 module score, Positive vs Negative (sanity):\n")
print(tapply(hyp_cells$HypoxiaSig1, hyp_cells$hypoxia, mean))
print(wilcox.test(HypoxiaSig1 ~ hypoxia, data = hyp_cells))

# --- A2. rank ALL genes by the hypoxia contrast (Positive vs Negative) ---
hyp <- subset(obj, cells = rownames(hyp_cells))
Idents(hyp) <- factor(hyp$hypoxia, levels = c("Negative","Positive"))
# full bidirectional ranking: no only.pos, no LFC threshold, keep broad coverage
deg <- FindMarkers(hyp, ident.1 = "Positive", ident.2 = "Negative",
                   only.pos = FALSE, logfc.threshold = 0, min.pct = 0.05)
deg$gene <- rownames(deg)
write.csv(deg, "hypoxia_DE_full.csv", row.names = FALSE)

# rank metric = avg_log2FC (transparent, monotonic in effect direction).
ranks <- deg$avg_log2FC
names(ranks) <- deg$gene
ranks <- sort(ranks[!is.na(ranks) & is.finite(ranks)], decreasing = TRUE)
cat(sprintf("\n[A2] Ranked %d genes for GSEA (top5: %s | bottom5: %s)\n",
            length(ranks),
            paste(head(names(ranks),5), collapse=","),
            paste(tail(names(ranks),5), collapse=",")))

# --- A3. Hallmark gene sets (msigdbr API changed across versions) ---
msig <- tryCatch(msigdbr(species = "Homo sapiens", collection = "H"),
                 error = function(e) msigdbr(species = "Homo sapiens", category = "H"))
hallmark <- split(msig$gene_symbol, msig$gs_name)
# add the curated HIF core as its own testable set
hallmark[["HIF_CORE_HANDOFF"]] <- hif_core

# --- A4. run fgsea ---
fg <- fgsea(pathways = hallmark, stats = ranks, minSize = 5, maxSize = 500, eps = 0)
fg <- fg[order(-fg$NES), ]
fg_out <- fg
fg_out$leadingEdge <- vapply(fg_out$leadingEdge, paste, collapse = ";", FUN.VALUE = "")
write.csv(fg_out, "gsea_hallmark_results.csv", row.names = FALSE)

cat("\n[A4] Top 10 enriched Hallmark sets (by NES):\n")
print(as.data.frame(fg[order(-fg$NES), .(pathway, NES, pval, padj, size)][1:10, ]))

cat("\n[A4] HALLMARK_HYPOXIA + curated HIF core:\n")
print(as.data.frame(fg[pathway %in% c("HALLMARK_HYPOXIA","HIF_CORE_HANDOFF"),
                       .(pathway, NES, pval, padj, size)]))

# --- A5. leading edge of HALLMARK_HYPOXIA = the tight, defensible core ---
le <- fg[pathway == "HALLMARK_HYPOXIA", ]$leadingEdge[[1]]
cat(sprintf("\n[A5] HALLMARK_HYPOXIA leading edge (%d genes) -- the tight core:\n", length(le)))
print(le)
cat("\n   handoff HIF core recovered in leading edge:\n")
print(setNames(hif_core %in% le, hif_core))
write.csv(data.frame(gene = le), "hif_core_leadingedge.csv", row.names = FALSE)

# enrichment plot (best-effort)
tryCatch({
  p <- plotEnrichment(hallmark[["HALLMARK_HYPOXIA"]], ranks) +
         ggtitle(sprintf("HALLMARK_HYPOXIA  (NES=%.2f, padj=%.3g)",
                         fg[pathway=="HALLMARK_HYPOXIA"]$NES,
                         fg[pathway=="HALLMARK_HYPOXIA"]$padj))
  ggsave("gsea_hypoxia_enrichment.png", p, width = 6, height = 4, dpi = 150)
  cat("\n   saved gsea_hypoxia_enrichment.png\n")
}, error = function(e) cat("   (enrichment plot skipped:", conditionMessage(e), ")\n"))

# ============================================================
#  PART B — DE-CONFOUND H2b   (revised to the genes that ACTUALLY exist)
#
#  Discovery (discover_immune.R) showed the Smart-seq2 CTC data is sparse:
#  of the 11 genes in the original Neutrophil1 'recruitment' score, only
#  S100A9 survives QC. Tumor-secreted recruiters (CSF3/CSF2/CXCL1/2/5/6/8)
#  are NOT detected at all -> the chemokine-recruitment hypothesis is
#  untestable here, and the H2b signal is a single neutrophil-intrinsic gene.
#  So Part B (1) proves the confound and (2) salvage-tests the trustworthy
#  H2a (hypoxia score) against a WBC-content covariate.
# ============================================================
cat("\n\n##############  PART B : DE-CONFOUND H2b  ##############\n")
g <- rownames(obj)

# --- B0. how much of the original recruitment panel even exists? ---
orig_neutro <- c("CXCL8","CXCL1","CXCL2","CXCL3","CXCL5","CXCL6","CSF3","CSF2",
                 "S100A8","S100A9","IL6")
present_orig <- intersect(orig_neutro, g)
cat(sprintf("\n[B0] Original Neutrophil1 panel present: %d/%d  -> {%s}\n",
            length(present_orig), length(orig_neutro), paste(present_orig, collapse=",")))
cat("     Tumor-secreted recruiters (CSF3/CSF2/CXCL1..6) NOT detected -> H2b is S100A9 alone.\n")
obj$S100A9 <- FetchData(obj, vars = "S100A9")[, 1]
cat(sprintf("[B0] cor(Neutrophil1 score, S100A9 expression) = %.3f  (the score IS S100A9)\n",
            cor(obj$Neutrophil1, obj$S100A9)))

# --- B1. WBC-content covariate from leukocyte/myeloid IDENTITY genes present ---
wbc_identity <- intersect(c("CD14","CD68","LYZ","CST3","TYMP","LCP1","FUT4"), g)
cat(sprintf("\n[B1] WBC-content covariate from %d myeloid/leukocyte genes: {%s}\n",
            length(wbc_identity), paste(wbc_identity, collapse = ",")))
obj <- AddModuleScore(obj, features = list(wbc_identity), name = "WBCcontent",
                      ctrl = 50, seed = 1)

# circadian cells with a WBC-attachment label
cc <- obj@meta.data[!is.na(obj$wbc_status), ]
cc$wbc_bin     <- as.integer(cc$wbc_status == "with_WBC")
cc$dataset_id  <- droplevels(factor(cc$dataset_id))
cat(sprintf("Cells used: %d  (with_WBC=%d, no_WBC=%d)\n",
            nrow(cc), sum(cc$wbc_bin), sum(1 - cc$wbc_bin)))

# --- B2. the confound is real & measurable ---
cat("\n[B2] WBC-content (contamination proxy), with_WBC vs no_WBC:\n")
cat(sprintf("  with_WBC=%.3f  no_WBC=%.3f  wilcox p=%s\n",
    mean(cc$WBCcontent1[cc$wbc_bin==1]), mean(cc$WBCcontent1[cc$wbc_bin==0]),
    format.pval(wilcox.test(cc$WBCcontent1 ~ cc$wbc_bin)$p.value, digits = 3)))

# --- B3. is the H2b 'recruitment' signal just contamination? ---
cat("\n[B3] H2b: lm(Neutrophil1 ~ wbc_status [+ WBCcontent])  -- does the effect survive?\n")
b_raw <- lm(Neutrophil1 ~ wbc_bin, data = cc)
b_adj <- lm(Neutrophil1 ~ wbc_bin + WBCcontent1, data = cc)
cat("  wbc_status coef  RAW         : ", paste(round(coef(summary(b_raw))["wbc_bin", ], 4), collapse="  "), "\n")
cat("  wbc_status coef  +WBCcontent : ", paste(round(coef(summary(b_adj))["wbc_bin", ], 4), collapse="  "), "\n")
cat("  (Estimate  StdErr  t  p)\n")

# --- B4. SALVAGE: is the trustworthy H2a (hypoxia score) robust to contamination? ---
cat("\n[B4] H2a: lm(HypoxiaSig ~ wbc_status [+ WBCcontent] [+ donor])  -- robustness:\n")
cat(sprintf("   raw wilcox HypoxiaSig by wbc_status: p=%s\n",
    format.pval(wilcox.test(HypoxiaSig1 ~ wbc_bin, data = cc)$p.value, digits = 3)))
h_raw  <- lm(HypoxiaSig1 ~ wbc_bin, data = cc)
h_adj  <- lm(HypoxiaSig1 ~ wbc_bin + WBCcontent1, data = cc)
h_adjd <- lm(HypoxiaSig1 ~ wbc_bin + WBCcontent1 + dataset_id, data = cc)
cat("  wbc_status coef  RAW         : ", paste(round(coef(summary(h_raw))["wbc_bin", ], 4), collapse="  "), "\n")
cat("  wbc_status coef  +WBCcontent : ", paste(round(coef(summary(h_adj))["wbc_bin", ], 4), collapse="  "), "\n")
cat("  wbc_status coef  +WBC+donor  : ", paste(round(coef(summary(h_adjd))["wbc_bin", ], 4), collapse="  "), "\n")
cat(sprintf("  cor(HypoxiaSig, WBCcontent) = %.3f  (overlap check; want LOW => no circularity)\n",
    cor(cc$HypoxiaSig1, cc$WBCcontent1)))

# --- B5. save scores ---
write.csv(cc[, c("dataset_id","wbc_status","HypoxiaSig1","Neutrophil1","S100A9","WBCcontent1")],
          "deconfound_scores.csv", row.names = TRUE)

cat("\n\nINTERPRETATION (read B3 vs B4):\n")
cat("  * H2b recruitment = S100A9 only (a neutrophil-INTRINSIC gene). If its wbc effect\n")
cat("    collapses once WBCcontent is controlled -> it is CONTAMINATION; drop the claim.\n")
cat("  * H2a hypoxia score: if its wbc effect SURVIVES WBCcontent (+donor) -> robust,\n")
cat("    i.e. WBC-attached CTCs are genuinely more hypoxic, not an artifact.\n")

cat("\nDONE. Outputs:\n",
    " - hypoxia_DE_full.csv          (full ranked DE list)\n",
    " - gsea_hallmark_results.csv    (fgsea Hallmark table)\n",
    " - hif_core_leadingedge.csv     (TIGHT hypoxia core <- use this, not the 847)\n",
    " - gsea_hypoxia_enrichment.png  (enrichment plot)\n",
    " - deconfound_scores.csv        (per-cell scores)\n")
