# What immune / chemokine / leukocyte genes actually survive in ctc_obj.rds?
suppressMessages(library(Seurat))
obj <- readRDS("ctc_obj.rds"); g <- rownames(obj)

show <- function(label, pattern_or_list, is_pattern = TRUE) {
  hits <- if (is_pattern) grep(pattern_or_list, g, value = TRUE) else intersect(pattern_or_list, g)
  cat(sprintf("\n%-22s (%d): %s\n", label, length(hits), paste(sort(hits), collapse = ", ")))
}

cat("Total genes:", length(g), "\n")
show("CXCL* chemokines",  "^CXCL")
show("CCL* chemokines",   "^CCL[0-9]")
show("CXCR/CCR receptors","^C[XC]CR")
show("S100A*",            "^S100A")
show("IL* cytokines",     "^IL[0-9]")
show("CSF*",              "^CSF")

# the EXACT panel used to build the original Neutrophil1 score
orig_neutro <- c("CXCL8","CXCL1","CXCL2","CXCL3","CXCL5","CXCL6","CSF3","CSF2",
                 "S100A8","S100A9","IL6")
show("ORIG Neutrophil1 panel", orig_neutro, is_pattern = FALSE)

# candidate leukocyte / neutrophil identity markers
leuko <- c("PTPRC","CD45","FCGR3B","CSF3R","ELANE","MPO","LCN2","LTF","DEFA1","DEFA3",
           "FUT4","ITGAM","SELL","CXCR2","FPR1","MNDA","CEACAM3","FCN1","LST1","AIF1",
           "TYROBP","FCER1G","CD53","LCP1","CD14","LYZ","CD68","NCF1","NCF2","FGL2",
           "CSF1R","CD163","S100A12","CAMP","RETN","CST3","TYMP","SRGN","CYBB")
show("leukocyte/myeloid markers", leuko, is_pattern = FALSE)

# tumor / epithelial identity (to confirm these ARE tumor cells)
epi <- c("EPCAM","KRT8","KRT18","KRT19","KRT7","CDH1","KRT5","ELF3","CLDN3","CLDN4",
         "KRT17","ESR1","ERBB2","KRT81")
show("epithelial/tumor markers", epi, is_pattern = FALSE)
