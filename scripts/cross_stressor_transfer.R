# ============================================================
# CROSS-STRESSOR TRANSFER  —  decisive experiment
# Thesis: a universal CTC stress program transfers between
#         hypoxia (GSE126669) and circadian-rest (GSE180097)
#
# Protocol: leave-one-stressor-out zero-shot transfer
#   train on stressor A, predict stressor B, and vice versa
# Controls: within-domain CV (signal exists?) + label
#           permutation null (is transfer real or chance?)
# Output : AUROCs, permutation p-values, shared gene signature
# ============================================================

library(Seurat)

# ---- deps (install once if missing) ----
for (p in c("glmnet", "pROC")) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}
library(glmnet)
library(pROC)
set.seed(1)

# ============================================================
# LOAD labeled harmony object
# (has $hypoxia, $timepoint, $stress_condition, $dataset_id)
# ============================================================
load("post_integration_stress_v2.RData")
obj <- seurat_harmony
obj <- JoinLayers(obj)

# ---- expression matrix: variable features, log1p data layer ----
feats <- VariableFeatures(obj)
X <- t(as.matrix(GetAssayData(obj, layer = "data")[feats, ]))  # cells x genes
X <- X[, apply(X, 2, var) > 0]                                  # drop constant genes

# ============================================================
# DEFINE THE TWO DOMAINS  (1 = Stressed, 0 = Baseline)
# Domain A = hypoxia       : Positive=1, Negative=0
# Domain B = circadian     : resting =1, active  =0
# ============================================================
isA <- !is.na(obj$hypoxia)   & obj$hypoxia   %in% c("Positive", "Negative")
isB <- !is.na(obj$timepoint) & obj$timepoint %in% c("resting",  "active")

yA <- ifelse(obj$hypoxia[isA]   == "Positive", 1L, 0L)
yB <- ifelse(obj$timepoint[isB] == "resting",  1L, 0L)
XA <- X[isA, , drop = FALSE]
XB <- X[isB, , drop = FALSE]

cat(sprintf("Domain A (hypoxia)  : %d cells | Stressed %d / Baseline %d\n",
            nrow(XA), sum(yA), sum(1 - yA)))
cat(sprintf("Domain B (circadian): %d cells | Stressed %d / Baseline %d\n\n",
            nrow(XB), sum(yB), sum(1 - yB)))

# ---- helper: fit L1-logistic (lambda by CV), return fitted cv.glmnet ----
fit_l1 <- function(Xtr, ytr) {
  nf <- min(10, sum(ytr == 1), sum(ytr == 0))      # safe folds for tiny n
  cv.glmnet(Xtr, ytr, family = "binomial",
            type.measure = "auc", alpha = 1, nfolds = max(3, nf))
}
auc_of <- function(model, Xte, yte) {
  p <- as.numeric(predict(model, Xte, s = "lambda.min", type = "response"))
  as.numeric(pROC::auc(pROC::roc(yte, p, quiet = TRUE)))
}

# ============================================================
# (0) SANITY: within-domain CV AUROC — is there ANY signal?
# ============================================================
cat("== Within-domain (does a stress signal even exist?) ==\n")
mA <- fit_l1(XA, yA); cat(sprintf("  Hypoxia   CV AUROC: %.3f\n", max(mA$cvm)))
mB <- fit_l1(XB, yB); cat(sprintf("  Circadian CV AUROC: %.3f\n\n", max(mB$cvm)))

# ============================================================
# (1) THE DECISIVE TEST: zero-shot cross-stressor transfer
# ============================================================
auc_A2B <- auc_of(mA, XB, yB)   # train hypoxia  -> predict circadian
auc_B2A <- auc_of(mB, XA, yA)   # train circadian-> predict hypoxia
cat("== Cross-stressor zero-shot transfer ==\n")
cat(sprintf("  Hypoxia  -> Circadian AUROC: %.3f\n", auc_A2B))
cat(sprintf("  Circadian-> Hypoxia   AUROC: %.3f\n\n", auc_B2A))

# ============================================================
# (2) PERMUTATION NULL: shuffle train labels, re-transfer
#     proves the transfer is real, not batch/chance artifact
# ============================================================
nperm <- 200
permA2B <- replicate(nperm, auc_of(fit_l1(XA, sample(yA)), XB, yB))
permB2A <- replicate(nperm, auc_of(fit_l1(XB, sample(yB)), XA, yA))
pA2B <- (1 + sum(permA2B >= auc_A2B)) / (1 + nperm)
pB2A <- (1 + sum(permB2A >= auc_B2A)) / (1 + nperm)
cat("== Permutation p-values (label-shuffle null) ==\n")
cat(sprintf("  Hypoxia  -> Circadian: p = %.4f (null mean %.3f)\n", pA2B, mean(permA2B)))
cat(sprintf("  Circadian-> Hypoxia  : p = %.4f (null mean %.3f)\n\n", pB2A, mean(permB2A)))

# ============================================================
# (3) UNIVERSAL STRESS SIGNATURE
#     genes non-zero AND sign-concordant in BOTH domains
# ============================================================
coef_of <- function(m) { b <- as.matrix(coef(m, s = "lambda.min"))[-1, 1]; b }
bA <- coef_of(mA); bB <- coef_of(mB)
g  <- intersect(names(bA), names(bB))
shared <- g[bA[g] != 0 & bB[g] != 0 & sign(bA[g]) == sign(bB[g])]
sig <- data.frame(gene = shared,
                  weight_hypoxia   = bA[shared],
                  weight_circadian = bB[shared],
                  direction = ifelse(bA[shared] > 0, "UP_in_stress", "DOWN_in_stress"))
sig <- sig[order(-abs(sig$weight_hypoxia + sig$weight_circadian)), ]
cat(sprintf("== Universal CTC stress signature: %d genes ==\n", nrow(sig)))
print(utils::head(sig, 30), row.names = FALSE)
write.csv(sig, "universal_stress_signature.csv", row.names = FALSE)
cat("\nSaved -> universal_stress_signature.csv\n")

# ============================================================
# INTERPRET
#  - Within-domain AUROC >> 0.5  => a stress signal exists.
#  - BOTH transfer AUROCs > 0.5 with p < 0.05 => the program
#    is SHARED across orthogonal stressors  => thesis holds.
#  - 'universal_stress_signature.csv' = headline gene list;
#    next run GSEA on it (expect HIF + general stress: UPR,
#    NF-kB, p53, dormancy/EMT beyond hypoxia-only genes).
# ============================================================
