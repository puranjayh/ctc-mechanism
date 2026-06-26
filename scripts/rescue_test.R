# ============================================================
# RESCUE TESTS  —  is the stress program shared at a level
# more robust than single-gene linear transfer?
#   Test 1: logFC concordance (do same genes move same way?)
#   Test 2: pathway-activity cross-stressor transfer
# No new compiled deps (uses Seurat + base R + glmnet already installed)
# ============================================================

library(Seurat)
set.seed(1)

load("post_integration_stress_v2.RData")
obj <- JoinLayers(seurat_harmony)

# ---- domains: A = hypoxia (Pos/Neg), B = circadian (resting/active) ----
A <- subset(obj, cells = colnames(obj)[!is.na(obj$hypoxia)   & obj$hypoxia   %in% c("Positive","Negative")])
B <- subset(obj, cells = colnames(obj)[!is.na(obj$timepoint) & obj$timepoint %in% c("resting","active")])
Idents(A) <- ifelse(A$hypoxia   == "Positive", "Stressed", "Baseline")
Idents(B) <- ifelse(B$timepoint == "resting",  "Stressed", "Baseline")

# ============================================================
# TEST 1 — SIGNATURE CONCORDANCE
# per-gene logFC (Stressed vs Baseline) in each domain, correlate
# ============================================================
deA <- FindMarkers(A, ident.1 = "Stressed", ident.2 = "Baseline",
                   logfc.threshold = 0, min.pct = 0, test.use = "wilcox")
deB <- FindMarkers(B, ident.1 = "Stressed", ident.2 = "Baseline",
                   logfc.threshold = 0, min.pct = 0, test.use = "wilcox")

g  <- intersect(rownames(deA), rownames(deB))
lA <- deA[g, "avg_log2FC"]; lB <- deB[g, "avg_log2FC"]
keep <- is.finite(lA) & is.finite(lB)
lA <- lA[keep]; lB <- lB[keep]; g <- g[keep]

pe <- cor.test(lA, lB, method = "pearson")
sp <- cor.test(lA, lB, method = "spearman")
cat("== TEST 1: logFC concordance (hypoxia vs circadian) ==\n")
cat(sprintf("  genes compared : %d\n", length(g)))
cat(sprintf("  Pearson  r = %.3f  (p = %.2e)\n", pe$estimate, pe$p.value))
cat(sprintf("  Spearman r = %.3f  (p = %.2e)\n", sp$estimate, sp$p.value))

# top concordant genes (moved strongly in SAME direction in both)
conc <- data.frame(gene = g, logFC_hyp = lA, logFC_circ = lB,
                   score = sign(lA) * sign(lB) * (abs(lA) + abs(lB)))
conc <- conc[conc$logFC_hyp * conc$logFC_circ > 0, ]
conc <- conc[order(-conc$score), ]
cat("\n  Top 20 concordant stress genes:\n")
print(utils::head(conc[, 1:3], 20), row.names = FALSE)
write.csv(conc, "concordant_stress_signature.csv", row.names = FALSE)

# ============================================================
# TEST 2 — PATHWAY-ACTIVITY CROSS-STRESSOR TRANSFER
# score cells with curated stress programs, transfer on scores
# ============================================================
pw <- list(
  HYPOXIA   = c("VEGFA","CA9","SLC2A1","LDHA","PGK1","HK2","BNIP3","PDK1","ALDOA","ENO1","NDRG1"),
  UPR       = c("HSPA5","DDIT3","ATF4","XBP1","HERPUD1","EDEM1","DNAJB9","ATF3"),
  P53       = c("CDKN1A","MDM2","BAX","GADD45A","BBC3","SESN1"),
  NFKB      = c("NFKBIA","TNFAIP3","CXCL8","IL6","CCL2","BIRC3","RELB"),
  ISR       = c("ATF4","DDIT3","TRIB3","ASNS","CHAC1","PPP1R15A","EIF4EBP1"),
  DORMANCY  = c("CDKN1B","CDKN2A","NR2F1","SOX9","THBS1")
)
present <- rownames(obj)
for (nm in names(pw)) {
  feat <- intersect(pw[[nm]], present)
  obj  <- AddModuleScore(obj, features = list(feat), name = nm, ctrl = min(50, length(feat)*5))
}
scol <- paste0(names(pw), "1")

isA <- !is.na(obj$hypoxia)   & obj$hypoxia   %in% c("Positive","Negative")
isB <- !is.na(obj$timepoint) & obj$timepoint %in% c("resting","active")
yA  <- ifelse(obj$hypoxia[isA]   == "Positive", 1L, 0L)
yB  <- ifelse(obj$timepoint[isB] == "resting",  1L, 0L)
SA  <- as.data.frame(obj@meta.data[isA, scol]); SB <- as.data.frame(obj@meta.data[isB, scol])
colnames(SA) <- colnames(SB) <- names(pw)

cat("\n== TEST 2: pathway-score cross-stressor transfer ==\n")
# direction check: does each pathway move the SAME way in both domains?
cat("  Pathway shift (mean Stressed - mean Baseline):\n")
for (nm in names(pw)) {
  dA <- mean(SA[yA==1, nm]) - mean(SA[yA==0, nm])
  dB <- mean(SB[yB==1, nm]) - mean(SB[yB==0, nm])
  cat(sprintf("    %-9s hypoxia %+.3f | circadian %+.3f | %s\n",
              nm, dA, dB, ifelse(sign(dA)==sign(dB), "CONCORDANT", "discordant")))
}

# logistic transfer on the 6 pathway scores
library(pROC)
fitG <- function(S, y) glm(y ~ ., data = cbind(S, y = y), family = "binomial")
au <- function(m, S, y) as.numeric(pROC::auc(pROC::roc(y, predict(m, S, type="response"), quiet=TRUE)))
mA <- fitG(SA, yA); mB <- fitG(SB, yB)
cat(sprintf("\n  Hypoxia  -> Circadian (pathway) AUROC: %.3f\n", au(mA, SB, yB)))
cat(sprintf("  Circadian-> Hypoxia   (pathway) AUROC: %.3f\n", au(mB, SA, yA)))

# ============================================================
# VERDICT
#  TEST 1 Pearson/Spearman r > 0 with small p  => shared program
#         exists at the directional/gene level (universal stress).
#  TEST 2 'CONCORDANT' pathways + transfer AUROC > 0.6 => shared
#         program is pathway-encoded -> strong, interpretable result.
#  If both null: stress is stressor-specific -> reframe + motivates
#         a domain-invariant (GNN/adversarial) method as the fix.
# ============================================================
