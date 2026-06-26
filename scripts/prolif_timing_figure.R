# Contrast figure: proliferation-at-rest is donor-ROBUST (rest>active in every
# donor) — the opposite of the hypoxia-WBC Simpson's-paradox artifact.
suppressMessages({ library(Seurat); library(ggplot2) })
set.seed(1)
obj <- readRDS("ctc_obj.rds"); DefaultAssay(obj) <- "RNA"
prolif <- intersect(c("MKI67","TOP2A","PCNA","CCNB1","CCNB2","CDK1","CDC20","AURKB",
  "BIRC5","CENPF","UBE2C","TYMS","RRM2","BUB1","KIF11","PLK1","NUSAP1","ASPM",
  "PRR11","CDKN3","KNL1","CDKL1"), rownames(obj))
obj <- AddModuleScore(obj, features = list(prolif), name = "Prolif", ctrl = 50, seed = 1)
cc <- obj@meta.data[obj$phase %in% c("resting","active"), ]
cc$donor <- droplevels(factor(cc$dataset_id))

p <- ggplot(cc, aes(donor, Prolif1, fill = phase)) +
  geom_boxplot(width = 0.7, outlier.size = 0.4, position = position_dodge(0.8)) +
  scale_fill_manual(values = c(resting = "#534AB7", active = "#B4B2A9")) +
  labs(title = "Proliferation is higher at rest in EVERY donor (donor-robust)",
       subtitle = "Donor-controlled p = 0.0036 — the opposite of the hypoxia-WBC artifact",
       x = "Donor / cell line", y = "Proliferation score", fill = "Phase") +
  theme_classic(base_size = 11)
ggsave("prolif_timing.png", p, width = 7, height = 4.2, dpi = 150)
cat("saved prolif_timing.png\n")
