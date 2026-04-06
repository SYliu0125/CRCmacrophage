# =============================================================================
# 02_processing.R
# Normalization, PCA, Harmony batch correction, UMAP, cluster verification
# =============================================================================

library(Seurat)
library(harmony)
library(ggplot2)
library(patchwork)
library(cowplot)
library(viridis)
library(RColorBrewer)

set.seed(42)
cat("=== Phase 2: Processing & Clustering Verification ===\n\n")

result_dir <- "results/02_processing"
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

# --- 1. Load QC'd Seurat object ---
cat("Loading Seurat object...\n")
seu <- readRDS("results/01_qc/seu_full_qc.rds")
cat(sprintf("  Loaded: %d genes x %d cells\n", nrow(seu), ncol(seu)))

# --- 2. Normalization ---
cat("Running SCTransform normalization...\n")
seu <- SCTransform(seu, verbose = TRUE, variable.features.n = 3000)

# --- 3. PCA ---
cat("Running PCA...\n")
seu <- RunPCA(seu, npcs = 50, verbose = FALSE)

# Elbow plot
p_elbow <- ElbowPlot(seu, ndims = 50) +
  geom_vline(xintercept = 30, linetype = "dashed", color = "red") +
  labs(title = "PCA Elbow Plot (red=30 PCs)") +
  theme_cowplot()

ggsave(file.path(result_dir, "pca_elbow.pdf"), p_elbow, width = 6, height = 4)
ggsave(file.path(result_dir, "pca_elbow.png"), p_elbow, width = 6, height = 4, dpi = 150)

# --- 4. Harmony batch correction ---
cat("Running Harmony batch correction (by batchID)...\n")
seu <- RunHarmony(seu, group.by.vars = "batchID", reduction = "pca",
                  assay.use = "SCT", dims.use = 1:30)

# --- 5. UMAP ---
cat("Running UMAP...\n")
seu <- RunUMAP(seu, reduction = "harmony", dims = 1:30, verbose = FALSE)

# --- 6. Verify published clusters on UMAP ---
cat("Generating UMAP plots...\n")

# Define colors for top-level cell types
celltype_colors <- c(
  "Epi"    = "#E41A1C", "TNKILC" = "#377EB8", "Myeloid" = "#4DAF4A",
  "Plasma" = "#984EA3", "B"      = "#FF7F00", "Strom"   = "#A65628",
  "Mast"   = "#F781BF"
)

# UMAP by clTopLevel
p1 <- DimPlot(seu, reduction = "umap", group.by = "clTopLevel",
              cols = celltype_colors, pt.size = 0.05, raster = TRUE) +
  labs(title = "UMAP - Major Cell Types") +
  theme_cowplot() +
  guides(color = guide_legend(override.aes = list(size = 3)))

# UMAP by cl295v11SubShort
p2 <- DimPlot(seu, reduction = "umap", group.by = "cl295v11SubShort",
              pt.size = 0.05, raster = TRUE, label = TRUE, label.size = 2) +
  labs(title = "UMAP - Subclusters") +
  theme_cowplot() + NoLegend()

# UMAP by specimen type
p3 <- DimPlot(seu, reduction = "umap", group.by = "SPECIMEN_TYPE",
              pt.size = 0.05, raster = TRUE) +
  labs(title = "UMAP - Specimen Type") +
  theme_cowplot()

# UMAP by MMR status
p4 <- DimPlot(seu, reduction = "umap", group.by = "MMRStatus",
              pt.size = 0.05, raster = TRUE) +
  labs(title = "UMAP - MMR Status") +
  theme_cowplot()

# UMAP by tissue site
p5 <- DimPlot(seu, reduction = "umap", group.by = "TissueSiteSimple",
              pt.size = 0.05, raster = TRUE) +
  labs(title = "UMAP - Tissue Site") +
  theme_cowplot()

# UMAP by batch
p6 <- DimPlot(seu, reduction = "umap", group.by = "batchID",
              pt.size = 0.05, raster = TRUE) +
  labs(title = "UMAP - Batch (post-Harmony)") +
  theme_cowplot() + NoLegend()

# Save individual plots
ggsave(file.path(result_dir, "umap_celltype.pdf"), p1, width = 8, height = 6)
ggsave(file.path(result_dir, "umap_celltype.png"), p1, width = 8, height = 6, dpi = 150)
ggsave(file.path(result_dir, "umap_subclusters.pdf"), p2, width = 10, height = 8)
ggsave(file.path(result_dir, "umap_subclusters.png"), p2, width = 10, height = 8, dpi = 150)
ggsave(file.path(result_dir, "umap_specimen.pdf"), p3, width = 8, height = 6)
ggsave(file.path(result_dir, "umap_mmr.pdf"), p4, width = 8, height = 6)
ggsave(file.path(result_dir, "umap_tissue.pdf"), p5, width = 8, height = 6)
ggsave(file.path(result_dir, "umap_batch.pdf"), p6, width = 8, height = 6)

# Combined panel
p_combined <- (p1 | p2) / (p3 | p4) +
  plot_annotation(title = "UMAP Overview - CRC scRNA-seq (Pelka et al. 2021)")

ggsave(file.path(result_dir, "umap_combined.pdf"), p_combined, width = 16, height = 12)
ggsave(file.path(result_dir, "umap_combined.png"), p_combined, width = 16, height = 12, dpi = 150)

# --- 7. Feature plots for key markers ---
cat("Generating marker feature plots...\n")
markers <- c("EPCAM", "CD3E", "CD14", "CD68", "CD19", "MS4A1",
             "JCHAIN", "COL1A1", "KIT", "TPSAB1")

p_feat <- FeaturePlot(seu, features = markers, ncol = 5,
                       pt.size = 0.05, raster = TRUE) &
  scale_color_viridis(option = "magma") &
  theme_cowplot() &
  theme(plot.title = element_text(size = 10))

ggsave(file.path(result_dir, "umap_markers.pdf"), p_feat, width = 20, height = 8)
ggsave(file.path(result_dir, "umap_markers.png"), p_feat, width = 20, height = 8, dpi = 150)

# --- 8. Save processed object ---
cat("\nSaving processed Seurat object...\n")
saveRDS(seu, file.path(result_dir, "seu_full_processed.rds"))

cat("\n Phase 2 complete!\n")
