# Run this ONCE on a fresh machine before anything else.
# Takes ~10-20 min depending on internet speed.

options(repos = c(CRAN = "https://cloud.r-project.org"))

# CRAN packages
cran_pkgs <- c("Seurat", "SeuratObject", "ggplot2", "patchwork", "dplyr", "msigdbr", "BiocManager")
for (p in cran_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    message("Installing ", p, "...")
    install.packages(p)
  } else {
    message(p, " already installed.")
  }
}

# Bioconductor packages
bioc_pkgs <- c("fgsea", "org.Hs.eg.db", "AnnotationDbi")
for (p in bioc_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    message("Installing ", p, " from Bioconductor...")
    BiocManager::install(p, ask = FALSE)
  } else {
    message(p, " already installed.")
  }
}

message("\nAll packages installed. Run gsea_deconfound.R next.")
