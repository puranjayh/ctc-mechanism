# Per-donor breakdown of H1b (cluster) and H2a (WBC) — is either effect
# real within donors, or just a between-donor (cell-line) artifact?
suppressMessages(library(Seurat))
obj <- readRDS("ctc_obj.rds")
cc <- obj@meta.data[!is.na(obj$wbc_status), ]
cc$dataset_id <- droplevels(factor(cc$dataset_id))

cat("=== CONFOUNDING STRUCTURE ===\n")
cat("\nwbc_status x donor:\n");   print(table(cc$wbc_status, cc$dataset_id))
cat("\nclumping  x donor:\n");    print(table(cc$clumping,  cc$dataset_id))
cat("\nsample_type x donor:\n");  print(table(cc$sample_type, cc$dataset_id))

cat("\n=== per-donor baseline HypoxiaSig1 (mean) ===\n")
print(round(tapply(cc$HypoxiaSig1, cc$dataset_id, mean), 3))

split_test <- function(s, var, a, b) {
  na <- sum(s[[var]] == a); nb <- sum(s[[var]] == b)
  if (na < 3 || nb < 3) return(sprintf("with=%d no=%d -- insufficient", na, nb))
  p <- wilcox.test(s$HypoxiaSig1 ~ s[[var]])$p.value
  sprintf("%s=%.3f (n=%d)  %s=%.3f (n=%d)  p=%.3g",
          a, mean(s$HypoxiaSig1[s[[var]]==a]), na,
          b, mean(s$HypoxiaSig1[s[[var]]==b]), nb, p)
}

cat("\n=== H2a per donor: HypoxiaSig1  with_WBC vs no_WBC ===\n")
for (d in levels(cc$dataset_id))
  cat(sprintf("  %-8s %s\n", d, split_test(cc[cc$dataset_id==d,], "wbc_status", "with_WBC", "no_WBC")))

cat("\n=== H1b per donor: HypoxiaSig1  cluster vs single ===\n")
for (d in levels(cc$dataset_id))
  cat(sprintf("  %-8s %s\n", d, split_test(cc[cc$dataset_id==d,], "clumping", "cluster", "single")))

# Cliff's delta effect sizes (rank-based, robust)
cliff <- function(x, y) { o <- outer(x, y, "-"); (sum(o>0) - sum(o<0)) / length(o) }
cat("\n=== overall effect sizes (Cliff's delta, |>0.33| ~ medium) ===\n")
cat(sprintf("  H2a hypoxia  with_WBC vs no_WBC : delta = %+.3f\n",
    cliff(cc$HypoxiaSig1[cc$wbc_status=="with_WBC"], cc$HypoxiaSig1[cc$wbc_status=="no_WBC"])))
cat(sprintf("  H1b hypoxia  cluster  vs single : delta = %+.3f\n",
    cliff(cc$HypoxiaSig1[cc$clumping=="cluster"],   cc$HypoxiaSig1[cc$clumping=="single"])))

# H1b with donor control (the analogue of the H2a collapse test)
cat("\n=== H1b robustness: lm(HypoxiaSig ~ clumping [+ donor]) ===\n")
c_raw <- lm(HypoxiaSig1 ~ I(clumping=="cluster"), data = cc)
c_don <- lm(HypoxiaSig1 ~ I(clumping=="cluster") + dataset_id, data = cc)
cat("  cluster coef  RAW    : ", paste(round(coef(summary(c_raw))[2, ], 4), collapse="  "), "\n")
cat("  cluster coef  +donor : ", paste(round(coef(summary(c_don))[2, ], 4), collapse="  "), "\n")
cat("  (Estimate  StdErr  t  p)\n")
