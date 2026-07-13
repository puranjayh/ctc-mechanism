# ============================================================
#  Publication-quality figures  ->  figures/*.tiff (300 DPI, LZW) + *.pdf (vector)
#  Regenerates the 3 key figures at submission resolution.
# ============================================================
suppressMessages({
  library(Seurat); library(ggplot2); library(patchwork)
  library(fgsea); library(msigdbr); library(dplyr)
})
set.seed(1)
obj <- readRDS("ctc_obj.rds"); DefaultAssay(obj) <- "RNA"
dir.create("figures", showWarnings = FALSE)

save_pub <- function(plot, stem, w, h) {
  f <- file.path("figures", stem)
  ggsave(paste0(f, ".tiff"), plot, width = w, height = h, dpi = 300)       # 300 DPI raster (TIFF)
  ggsave(paste0(f, ".pdf"),  plot, width = w, height = h, device = "pdf")  # vector (base pdf; cairo unavailable)
  stopifnot(file.exists(paste0(f, ".tiff")), file.exists(paste0(f, ".pdf")))
  cat(sprintf("  %-24s  tiff + pdf OK\n", stem))
}

# ---- Fig 1. Simpson's paradox: hypoxia vs WBC attachment ----
cc <- obj@meta.data[!is.na(obj$wbc_status), ]; cc$donor <- droplevels(factor(cc$dataset_id))
pA <- ggplot(cc, aes(wbc_status, HypoxiaSig1, fill = wbc_status)) +
  geom_boxplot(width = 0.6, outlier.size = 0.4) +
  labs(title = "Pooled (the illusion)", subtitle = "with_WBC looks more hypoxic (p = 0.008)",
       x = NULL, y = "Hypoxia score") +
  theme_classic(base_size = 12) + theme(legend.position = "none")
pB <- ggplot(cc, aes(donor, HypoxiaSig1, fill = wbc_status)) +
  geom_boxplot(width = 0.7, outlier.size = 0.4, position = position_dodge(0.8)) +
  labs(title = "Within donor (the truth)", subtitle = "No within-donor gap (p = 0.62 / 0.86 / 0.84)",
       x = "Donor / cell line", y = "Hypoxia score", fill = "WBC") +
  theme_classic(base_size = 12)
save_pub(pA | pB, "fig1_simpson_hypoxia_wbc", 10, 4)

# ---- Fig 2. Proliferation is higher at rest, in every donor ----
prolif <- intersect(c("MKI67","TOP2A","PCNA","CCNB1","CCNB2","CDK1","CDC20","AURKB","BIRC5",
  "CENPF","UBE2C","TYMS","RRM2","BUB1","KIF11","PLK1","NUSAP1","ASPM","PRR11","CDKN3","KNL1","CDKL1"),
  rownames(obj))
obj <- AddModuleScore(obj, features = list(prolif), name = "Prolif", ctrl = 50, seed = 1)
cc2 <- obj@meta.data[obj$phase %in% c("resting","active"), ]; cc2$donor <- droplevels(factor(cc2$dataset_id))
p2 <- ggplot(cc2, aes(donor, Prolif1, fill = phase)) +
  geom_boxplot(width = 0.7, outlier.size = 0.4, position = position_dodge(0.8)) +
  scale_fill_manual(values = c(resting = "#534AB7", active = "#B4B2A9")) +
  labs(title = "Proliferation is higher at rest in every donor",
       subtitle = "Donor-controlled p = 0.0036", x = "Donor / cell line",
       y = "Proliferation score", fill = "Phase") +
  theme_classic(base_size = 12)
save_pub(p2, "fig2_proliferation_rest", 7, 4.2)

# ---- Fig 3. GSEA: HALLMARK_HYPOXIA enrichment ----
hyp <- subset(obj, cells = rownames(obj@meta.data[obj$hypoxia %in% c("Positive","Negative"), ]))
Idents(hyp) <- factor(hyp$hypoxia, levels = c("Negative","Positive"))
deg <- FindMarkers(hyp, ident.1 = "Positive", ident.2 = "Negative",
                   only.pos = FALSE, logfc.threshold = 0, min.pct = 0.05)
ranks <- sort(setNames(deg$avg_log2FC, rownames(deg)), decreasing = TRUE); ranks <- ranks[is.finite(ranks)]
msig <- tryCatch(msigdbr(species = "Homo sapiens", collection = "H"),
                 error = function(e) msigdbr(species = "Homo sapiens", category = "H"))
hallmark <- split(msig$gene_symbol, msig$gs_name)
fg <- fgsea(hallmark, ranks, minSize = 5, maxSize = 500, eps = 0)
nes <- fg[fg$pathway == "HALLMARK_HYPOXIA", ]$NES; padj <- fg[fg$pathway == "HALLMARK_HYPOXIA", ]$padj
p3 <- plotEnrichment(hallmark[["HALLMARK_HYPOXIA"]], ranks) +
  labs(title = "HALLMARK_HYPOXIA enrichment",
       subtitle = sprintf("NES = %.2f, adjusted p = %.3g", nes, padj),
       x = "Gene rank (hypoxia Positive vs Negative)", y = "Enrichment score")
save_pub(p3, "fig3_gsea_hypoxia", 6, 4)

cat("\nDONE. Publication figures in ./figures/  (300 DPI TIFF + vector PDF)\n")
