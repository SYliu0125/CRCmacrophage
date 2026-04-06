# =============================================================================
# 04_myeloid_analysis.R
# Deep-dive into myeloid/macrophage cells: subset, re-cluster, markers
# =============================================================================

library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)
library(cowplot)
library(viridis)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)

set.seed(42)
cat("=== Phase 4: Myeloid / Macrophage Deep-Dive ===\n\n")

result_dir <- "results/04_myeloid"
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

# --- 1. Load and subset ---
cat("Loading processed Seurat object...\n")
seu_full <- readRDS("results/02_processing/seu_full_processed.rds")

cat("Subsetting myeloid cells...\n")
seu <- subset(seu_full, subset = clTopLevel == "Myeloid")
cat(sprintf("  Myeloid subset: %d genes x %d cells\n", nrow(seu), ncol(seu)))
rm(seu_full); gc()

# --- 2. Re-process myeloid cells ---
cat("Re-processing myeloid cells...\n")
DefaultAssay(seu) <- "RNA"
seu <- NormalizeData(seu, verbose = FALSE)
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 3000, verbose = FALSE)
seu <- ScaleData(seu, verbose = FALSE)
seu <- RunPCA(seu, npcs = 30, verbose = FALSE)

# Harmony on myeloid subset
library(harmony)
seu <- RunHarmony(seu, group.by.vars = "batchID", reduction = "pca", dims.use = 1:30)
seu <- RunUMAP(seu, reduction = "harmony", dims = 1:20, verbose = FALSE)

# --- 3. Myeloid subcluster UMAP ---
cat("Generating myeloid UMAP plots...\n")

myeloid_colors <- c(
  "cM01 (Monocyte)"       = "#1f77b4",
  "cM02 (Macrophage-like)" = "#d62728",
  "cM03 (DC1)"             = "#2ca02c",
  "cM04 (DC2)"             = "#ff7f0e",
  "cM05 (DC2 C1Q+)"        = "#9467bd",
  "cM06 (DC IL22RA2)"      = "#8c564b",
  "cM07 (pDC)"             = "#e377c2",
  "cM08 (AS-DC)"           = "#7f7f7f",
  "cM09 (mregDC)"          = "#bcbd22",
  "cM10 (Granulocyte)"     = "#17becf"
)

p1 <- DimPlot(seu, reduction = "umap", group.by = "cl295v11SubFull",
              cols = myeloid_colors, pt.size = 0.3, label = TRUE, label.size = 3) +
  labs(title = "Myeloid UMAP - Subclusters") +
  theme_cowplot()

ggsave(file.path(result_dir, "myeloid_umap_subclusters.pdf"), p1, width = 10, height = 7)
ggsave(file.path(result_dir, "myeloid_umap_subclusters.png"), p1, width = 10, height = 7, dpi = 200)

# UMAP by specimen type
p2 <- DimPlot(seu, reduction = "umap", group.by = "SPECIMEN_TYPE",
              pt.size = 0.2) +
  labs(title = "Myeloid UMAP - Specimen Type") +
  theme_cowplot()

# UMAP by MMR status
p3 <- DimPlot(seu, reduction = "umap", group.by = "MMRStatus",
              pt.size = 0.2) +
  labs(title = "Myeloid UMAP - MMR Status") +
  theme_cowplot()

p_combined <- p1 / (p2 | p3)
ggsave(file.path(result_dir, "myeloid_umap_combined.pdf"), p_combined, width = 12, height = 14)
ggsave(file.path(result_dir, "myeloid_umap_combined.png"), p_combined, width = 12, height = 14, dpi = 150)

# --- 4. Key macrophage marker expression ---
cat("Plotting macrophage markers...\n")

# Key markers
macro_markers <- c(
  # Pan-myeloid
  "CD14", "CD68", "CSF1R", "FCGR3A",
  # M1/pro-inflammatory
  "TNF", "IL1B", "CXCL9", "CXCL10", "IDO1", "NOS2",
  # M2/anti-inflammatory
  "CD163", "MRC1", "MSR1", "TGFBI", "IL10",
  # SPP1+ macrophages
  "SPP1", "VEGFA", "MARCO",
  # C1Q+ macrophages
  "C1QA", "C1QB", "C1QC", "APOE",
  # DC markers
  "CLEC9A", "XCR1", "CD1C", "FCER1A", "IRF7", "LILRA4",
  # Monocyte
  "S100A8", "S100A9", "VCAN", "FCN1",
  # Granulocyte
  "CSF3R", "CXCR2"
)

# Filter to genes that exist in the dataset
macro_markers <- macro_markers[macro_markers %in% rownames(seu)]

# Dot plot
p_dot <- DotPlot(seu, features = macro_markers, group.by = "cl295v11SubFull",
                 dot.scale = 6) +
  coord_flip() +
  scale_color_viridis(option = "magma") +
  labs(title = "Macrophage / Myeloid Marker Expression") +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.text.y = element_text(size = 8))

ggsave(file.path(result_dir, "myeloid_dotplot_markers.pdf"), p_dot, width = 10, height = 12)
ggsave(file.path(result_dir, "myeloid_dotplot_markers.png"), p_dot, width = 10, height = 12, dpi = 200)

# Feature plots for key markers
key_markers <- c("CD14", "CD68", "SPP1", "C1QA", "S100A8", "CD163",
                 "TNF", "IL1B", "VEGFA", "CLEC9A", "IRF7", "CSF3R")
key_markers <- key_markers[key_markers %in% rownames(seu)]

p_feat <- FeaturePlot(seu, features = key_markers, ncol = 4,
                       pt.size = 0.1, raster = TRUE) &
  scale_color_viridis(option = "magma") &
  theme_cowplot() &
  theme(plot.title = element_text(size = 10))

ggsave(file.path(result_dir, "myeloid_featureplots.pdf"), p_feat, width = 16, height = 12)
ggsave(file.path(result_dir, "myeloid_featureplots.png"), p_feat, width = 16, height = 12, dpi = 150)

# --- 5. Find markers for each subcluster ---
cat("Finding markers for each myeloid subcluster...\n")
Idents(seu) <- "cl295v11SubFull"
myeloid_markers <- FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5,
                                   test.use = "wilcox", verbose = TRUE)

write.csv(myeloid_markers, file.path(result_dir, "myeloid_subcluster_markers.csv"), row.names = FALSE)

# Top 10 markers per cluster heatmap
top10 <- myeloid_markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC) %>%
  ungroup()

# Average expression heatmap
avg_expr <- AverageExpression(seu, features = unique(top10$gene), group.by = "cl295v11SubFull")
avg_mat <- as.matrix(avg_expr$RNA)
avg_mat_scaled <- t(scale(t(avg_mat)))

pdf(file.path(result_dir, "myeloid_top_markers_heatmap.pdf"), width = 10, height = 14)
ht <- Heatmap(avg_mat_scaled, name = "Z-score",
              col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
              row_names_gp = gpar(fontsize = 6),
              column_names_gp = gpar(fontsize = 8),
              cluster_columns = TRUE, cluster_rows = TRUE,
              column_title = "Top Markers per Myeloid Subcluster")
draw(ht)
dev.off()

# --- 6. M1/M2 scoring ---
cat("Computing M1/M2 macrophage scores...\n")

m1_genes <- c("TNF", "IL1B", "IL6", "CXCL9", "CXCL10", "CXCL11",
              "IDO1", "NOS2", "CD80", "CD86", "HLA-DRA")
m2_genes <- c("CD163", "MRC1", "MSR1", "TGFBI", "IL10", "ARG1",
              "STAB1", "FOLR2", "LYVE1", "CCL18")

m1_genes <- m1_genes[m1_genes %in% rownames(seu)]
m2_genes <- m2_genes[m2_genes %in% rownames(seu)]

seu <- AddModuleScore(seu, features = list(m1_genes), name = "M1_score")
seu <- AddModuleScore(seu, features = list(m2_genes), name = "M2_score")

# M1/M2 score plots
p_m1 <- FeaturePlot(seu, features = "M1_score1", pt.size = 0.2, raster = TRUE) +
  scale_color_viridis(option = "inferno") +
  labs(title = "M1 (Pro-inflammatory) Score") + theme_cowplot()

p_m2 <- FeaturePlot(seu, features = "M2_score1", pt.size = 0.2, raster = TRUE) +
  scale_color_viridis(option = "viridis") +
  labs(title = "M2 (Anti-inflammatory) Score") + theme_cowplot()

p_m1m2 <- p_m1 | p_m2
ggsave(file.path(result_dir, "m1_m2_scores_umap.pdf"), p_m1m2, width = 14, height = 6)
ggsave(file.path(result_dir, "m1_m2_scores_umap.png"), p_m1m2, width = 14, height = 6, dpi = 150)

# M1 vs M2 violin by subcluster
p_vln <- VlnPlot(seu, features = c("M1_score1", "M2_score1"),
                  group.by = "cl295v11SubFull", pt.size = 0, ncol = 2) &
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))

ggsave(file.path(result_dir, "m1_m2_violin_by_subcluster.pdf"), p_vln, width = 14, height = 5)
ggsave(file.path(result_dir, "m1_m2_violin_by_subcluster.png"), p_vln, width = 14, height = 5, dpi = 150)

# --- 7. Save myeloid Seurat object ---
cat("\nSaving myeloid Seurat object...\n")
saveRDS(seu, file.path(result_dir, "seu_myeloid.rds"))

cat("\n Phase 4 complete!\n")
