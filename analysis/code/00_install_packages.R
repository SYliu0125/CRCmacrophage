# =============================================================================
# 00_install_packages.R
# Install all required packages for CRC single-cell reanalysis
# Includes: Seurat, trajectory (monocle3, slingshot), SCENIC, and more
# =============================================================================

cat("=== Installing required packages for CRC scRNA-seq reanalysis ===\n\n")

# --- CRAN packages ---
cran_pkgs <- c(
  "Seurat", "SeuratObject", "hdf5r",
  "harmony", "patchwork", "cowplot", "gridExtra",
  "reshape2", "viridis", "igraph", "uwot", "Rtsne",
  "tidyverse", "data.table", "scales", "RColorBrewer",
  "pheatmap", "devtools", "R.utils"
)

missing_cran <- cran_pkgs[!sapply(cran_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_cran) > 0) {
  cat("Installing CRAN packages:", paste(missing_cran, collapse = ", "), "\n")
  install.packages(missing_cran, repos = "https://cloud.r-project.org", dependencies = TRUE)
} else {
  cat("All CRAN packages already installed.\n")
}

# --- Bioconductor packages ---
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

bioc_pkgs <- c(
  "SingleCellExperiment", "scater", "scran",
  "MAST", "clusterProfiler", "org.Hs.eg.db", "enrichplot",
  "ComplexHeatmap", "EnhancedVolcano", "msigdbr",
  "slingshot", "tradeSeq",
  "AUCell", "RcisTarget", "GENIE3",
  "SCENIC"
)

missing_bioc <- bioc_pkgs[!sapply(bioc_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_bioc) > 0) {
  cat("Installing Bioconductor packages:", paste(missing_bioc, collapse = ", "), "\n")
  BiocManager::install(missing_bioc, update = FALSE, ask = FALSE)
}

# --- GitHub packages (monocle3) ---
if (!requireNamespace("monocle3", quietly = TRUE)) {
  cat("Installing monocle3 from GitHub...\n")
  devtools::install_github("cole-trapnell-lab/monocle3", upgrade = "never")
}

# If SCENIC not available via Bioconductor, try GitHub
if (!requireNamespace("SCENIC", quietly = TRUE)) {
  cat("Installing SCENIC from GitHub...\n")
  devtools::install_github("aertslab/SCENIC", upgrade = "never")
}

# --- Verify all packages ---
cat("\n=== Verification ===\n")
all_pkgs <- c(cran_pkgs, bioc_pkgs, "monocle3")
status <- sapply(all_pkgs, requireNamespace, quietly = TRUE)
for (i in seq_along(all_pkgs)) {
  cat(sprintf("  %-25s %s\n", all_pkgs[i], ifelse(status[i], "OK", "FAILED")))
}

if (all(status)) {
  cat("\n All packages installed successfully!\n")
} else {
  failed <- all_pkgs[!status]
  cat("\n Some packages failed:", paste(failed, collapse = ", "), "\n")
  cat("  You may need to install these manually.\n")
}
