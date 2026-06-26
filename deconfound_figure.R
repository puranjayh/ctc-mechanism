# Simpson's-paradox figure: the pooled "hypoxia is higher in WBC-attached CTCs"
# signal is a between-donor artifact.
suppressMessages({ library(Seurat); library(ggplot2); library(patchwork) })
obj <- readRDS("ctc_obj.rds")
cc <- obj@meta.data[!is.na(obj$wbc_status), ]
cc$donor <- droplevels(factor(cc$dataset_id))

pA <- ggplot(cc, aes(wbc_status, HypoxiaSig1, fill = wbc_status)) +
  geom_boxplot(width = 0.6, outlier.size = 0.4) +
  labs(title = "Pooled (the illusion)",
       subtitle = "with_WBC looks more hypoxic  (wilcox p = 0.008)",
       x = NULL, y = "Hypoxia score") +
  theme_classic(base_size = 11) + theme(legend.position = "none")

pB <- ggplot(cc, aes(donor, HypoxiaSig1, fill = wbc_status)) +
  geom_boxplot(width = 0.7, outlier.size = 0.4, position = position_dodge(0.8)) +
  labs(title = "Within donor (the truth)",
       subtitle = "No within-donor gap; donors differ in baseline (p = 0.62 / 0.86 / 0.84)",
       x = "Donor / cell line", y = "Hypoxia score", fill = "WBC") +
  theme_classic(base_size = 11)

ggsave("deconfound_simpson.png", pA | pB, width = 10, height = 4, dpi = 150)
cat("saved deconfound_simpson.png\n")
