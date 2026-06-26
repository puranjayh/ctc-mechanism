# ============================================================
#  METASTASIS-TIMING (advisor request: "check in context of metastasis")
#  GSE180097 = Aceto sleep-metastasis data; rest phase is the pro-metastatic
#  window. Do the METASTATIC PHENOTYPES (clustering, proliferation) peak at
#  rest once DONOR is controlled?  And is hypoxia involved (it shouldn't be).
# ============================================================
suppressMessages({ library(Seurat); library(ggplot2) })
set.seed(1)
obj <- readRDS("ctc_obj.rds")
DefaultAssay(obj) <- "RNA"

cc <- obj@meta.data[obj$phase %in% c("resting","active"), ]
cc$donor <- droplevels(factor(cc$dataset_id))
cc$rest  <- as.integer(cc$phase == "resting")
cat("Cells:", nrow(cc), " | phase x donor:\n"); print(table(cc$phase, cc$donor))

# proliferation score (cell-cycle / mitotic) -----------------
prolif <- c("MKI67","TOP2A","PCNA","CCNB1","CCNB2","CDK1","CDC20","AURKB","BIRC5",
            "CENPF","UBE2C","TYMS","RRM2","BUB1","KIF11","PLK1","NUSAP1","ASPM",
            "PRR11","CDKN3","KNL1","CDKL1")
prolif <- intersect(prolif, rownames(obj))
cat(sprintf("\nProliferation panel present: %d genes {%s}\n", length(prolif), paste(prolif, collapse=",")))
obj <- AddModuleScore(obj, features = list(prolif), name = "Prolif", ctrl = 50, seed = 1)
cc$Prolif1 <- obj$Prolif1[rownames(cc)]

donor_test <- function(var) {
  cat(sprintf("\n--- %s : resting vs active ---\n", var))
  cat(sprintf("  pooled wilcox p = %s  (rest=%.3f, active=%.3f)\n",
      format.pval(wilcox.test(cc[[var]] ~ cc$rest)$p.value, digits=3),
      mean(cc[[var]][cc$rest==1]), mean(cc[[var]][cc$rest==0])))
  m1 <- lm(cc[[var]] ~ cc$rest); m2 <- lm(cc[[var]] ~ cc$rest + cc$donor)
  cat("  rest coef  RAW   :", paste(round(coef(summary(m1))[2,],4), collapse="  "), "\n")
  cat("  rest coef  +donor:", paste(round(coef(summary(m2))[2,],4), collapse="  "), "\n")
  for (d in levels(cc$donor)) { s <- cc[cc$donor==d,]
    if (length(unique(s$rest))==2 && min(table(s$rest))>=3)
      cat(sprintf("    %-8s rest=%.3f active=%.3f  p=%.3g\n", d,
          mean(s[[var]][s$rest==1]), mean(s[[var]][s$rest==0]),
          wilcox.test(s[[var]] ~ s$rest)$p.value)) }
}

# 1. PROLIFERATION timing (expect rest > active per Aceto)
donor_test("Prolif1")
# 2. HYPOXIA timing (expect flat; advisor's bridge)
donor_test("HypoxiaSig1")

# 3. CLUSTERING frequency: are CTCs more clustered at rest? (donor-controlled)
cat("\n--- clustering frequency : resting vs active ---\n")
print(table(cc$phase, cc$clumping))
gl <- glm(I(clumping=="cluster") ~ rest + donor, data = cc, family = binomial)
cat("  logistic rest coef (log-odds of clustering at rest), +donor:\n")
print(round(coef(summary(gl))["rest",,drop=FALSE], 4))
for (d in levels(cc$donor)) { s <- cc[cc$donor==d,]
  tb <- table(s$rest, s$clumping=="cluster")
  if (all(dim(tb)==c(2,2)))
    cat(sprintf("    %-8s fisher p=%.3g  (rest %%cluster=%.0f%%, active=%.0f%%)\n", d,
        fisher.test(tb)$p.value,
        100*mean(s$clumping[s$rest==1]=="cluster"), 100*mean(s$clumping[s$rest==0]=="cluster"))) }

cat("\nDONE. Read: does rest drive proliferation/clustering AFTER donor control?\n")
