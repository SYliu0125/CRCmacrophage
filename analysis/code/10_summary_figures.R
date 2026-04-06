# =============================================================================
# 10_summary_figures.R
# Publication-quality summary figures combining key results
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
library(gridExtra)

set.seed(42)
cat("=== Phase 10: Summary Figures ===\n\n")

result_dir <- "results/10_summary"
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

# --- Load objects ---
cat("Loading Seurat objects...\n")
seu_full    <- readRDS("results/02_processing/seu_full_processed.rds")
seu_myeloid <- readRDS("results/04_myeloid/seu_myeloid.rds")

# Colors
celltype_colors <- c(
  "Epi" = "#E41A1C", "TNKILC" = "#377EB8", "Myeloid" = "#4DAF4A",
  "Plasma" = "#984EA3", "B" = "#FF7F00", "Strom" = "#A65628", "Mast" = "#F781BF"
)

myeloid_colors <- c(
  "cM01 (Monocyte)" = "#1f77b4", "cM02 (Macrophage-like)" = "#d62728",
  "cM03 (DC1)" = "#2ca02c", "cM04 (DC2)" = "#ff7f0e",
  "cM05 (DC2 C1Q+)" = "#9467bd", "cM06 (DC IL22RA2)" = "#8c564b",
  "cM07 (pDC)" = "#e377c2", "cM08 (AS-DC)" = "#7f7f7f",
  "cM09 (mregDC)" = "#bcbd22", "cM10 (Granulocyte)" = "#17becf"
)

# =========================================================================
# FIGURE 1: Overview (Full dataset UMAP + cell type bar + markers)
# =========================================================================
cat("Generating Figure 1: Dataset Overview...\n")

p1a <- DimPlot(seu_full, reduction = "umap", group.by = "clTopLevel",
               cols = celltype_colors, pt.size = 0.02, raster = TRUE) +
  labs(title = "A. UMAP - All Cell Types") +
  theme_cowplot() +
  guides(color = guide_legend(override.aes = list(size = 3)))

# Cell counts bar
ct_df <- data.frame(table(seu_full$clTopLevel)) %>%
  arrange(desc(Freq)) %>%
  mutate(Var1 = factor(Var1, levels = Var1))

p1b <- ggplot(ct_df, aes(x = Var1, y = Freq / 1000, fill = Var1)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = celltype_colors) +
  geom_text(aes(label = scales::comma(Freq)), vjust = -0.3, size = 2.8) +
  labs(title = "B. Cell Type Counts", x = "", y = "Cells (×1000)") +
  theme_cowplot() + NoLegend() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 9))

# Key markers
key_genes <- c("EPCAM", "CD3E", "CD14", "JCHAIN", "MS4A1", "COL1A1")
key_genes <- key_genes[key_genes %in% rownames(seu_full)]

p1c <- FeaturePlot(seu_full, features = key_genes, ncol = 3,
                    pt.size = 0.02, raster = TRUE) &
  scale_color_viridis(option = "magma") &
  theme_cowplot() &
  theme(plot.title = element_text(size = 10),
        legend.key.size = unit(0.3, "cm"))

fig1 <- (p1a | p1b) / p1c +
  plot_annotation(title = "Figure 1: CRC scRNA-seq Dataset Overview (Pelka et al. 2021 Reanalysis)",
                  theme = theme(plot.title = element_text(size = 14, face = "bold")))

ggsave(file.path(result_dir, "figure1_overview.pdf"), fig1, width = 16, height = 14)
ggsave(file.path(result_dir, "figure1_overview.png"), fig1, width = 16, height = 14, dpi = 200)

rm(seu_full); gc()

# =========================================================================
# FIGURE 2: Myeloid Deep-Dive
# =========================================================================
cat("Generating Figure 2: Myeloid Deep-Dive...\n")

p2a <- DimPlot(seu_myeloid, reduction = "umap", group.by = "cl295v11SubFull",
               cols = myeloid_colors, pt.size = 0.3, label = TRUE, label.size = 3) +
  labs(title = "A. Myeloid UMAP - Subclusters") +
  theme_cowplot() + NoLegend()

# M1/M2 scores
if ("M1_score1" %in% colnames(seu_myeloid@meta.data)) {
  p2b <- FeaturePlot(seu_myeloid, features = "M1_score1", pt.size = 0.2, raster = TRUE) +
    scale_color_viridis(option = "inferno") +
    labs(title = "B. M1 Score") + theme_cowplot()
  
  p2c <- FeaturePlot(seu_myeloid, features = "M2_score1", pt.size = 0.2, raster = TRUE) +
    scale_color_viridis(option = "viridis") +
    labs(title = "C. M2 Score") + theme_cowplot()
} else {
  # Compute if not available
  m1_genes <- c("TNF", "IL1B", "IL6", "CXCL9", "CXCL10", "IDO1")
  m2_genes <- c("CD163", "MRC1", "MSR1", "TGFBI", "IL10", "STAB1")
  m1_genes <- m1_genes[m1_genes %in% rownames(seu_myeloid)]
  m2_genes <- m2_genes[m2_genes %in% rownames(seu_myeloid)]
  
  seu_myeloid <- AddModuleScore(seu_myeloid, features = list(m1_genes), name = "M1_score")
  seu_myeloid <- AddModuleScore(seu_myeloid, features = list(m2_genes), name = "M2_score")
  
  p2b <- FeaturePlot(seu_myeloid, features = "M1_score1", pt.size = 0.2, raster = TRUE) +
    scale_color_viridis(option = "inferno") +
    labs(title = "B. M1 Score") + theme_cowplot()
  
  p2c <- FeaturePlot(seu_myeloid, features = "M2_score1", pt.size = 0.2, raster = TRUE) +
    scale_color_viridis(option = "viridis") +
    labs(title = "C. M2 Score") + theme_cowplot()
}

# Key macrophage markers
mac_genes <- c("CD14", "CD68", "SPP1", "C1QA", "CD163", "S100A8")
mac_genes <- mac_genes[mac_genes %in% rownames(seu_myeloid)]

p2d <- FeaturePlot(seu_myeloid, features = mac_genes, ncol = 3,
                    pt.size = 0.1, raster = TRUE) &
  scale_color_viridis(option = "magma") &
  theme_cowplot() &
  theme(plot.title = element_text(size = 10))

fig2 <- (p2a | p2b | p2c) / p2d +
  plot_annotation(title = "Figure 2: Myeloid / Macrophage Analysis",
                  theme = theme(plot.title = element_text(size = 14, face = "bold")))

ggsave(file.path(result_dir, "figure2_myeloid.pdf"), fig2, width = 18, height = 14)
ggsave(file.path(result_dir, "figure2_myeloid.png"), fig2, width = 18, height = 14, dpi = 200)

# =========================================================================
# FIGURE 3: Trajectory & Regulon Summary
# =========================================================================
cat("Generating Figure 3: Trajectory & Regulon...\n")

# Try loading trajectory/regulon results
fig3_plots <- list()

# Trajectory pseudotime
if ("slingshot_pt1" %in% colnames(seu_myeloid@meta.data)) {
  fig3_plots[["trajectory"]] <- FeaturePlot(seu_myeloid, features = "slingshot_pt1",
                                             pt.size = 0.2, raster = TRUE) +
    scale_color_viridis(option = "inferno", na.value = "grey80") +
    labs(title = "A. Slingshot Pseudotime (Lineage 1)") + theme_cowplot()
}

# Regulon scores
reg_features <- grep("^Reg_", colnames(seu_myeloid@meta.data), value = TRUE)
if (length(reg_features) >= 4) {
  plot_regs <- head(reg_features, 4)
  for (i in seq_along(plot_regs)) {
    reg_name <- gsub("^Reg_(.+)1$", "\\1", plot_regs[i])
    reg_name <- gsub("_", " ", reg_name)
    fig3_plots[[paste0("reg_", i)]] <- FeaturePlot(seu_myeloid, features = plot_regs[i],
                                                     pt.size = 0.2, raster = TRUE) +
      scale_color_viridis(option = "magma") +
      labs(title = paste0(LETTERS[i + 1], ". ", reg_name, " regulon")) + theme_cowplot()
  }
}

if (length(fig3_plots) > 0) {
  fig3 <- wrap_plots(fig3_plots, ncol = 3) +
    plot_annotation(title = "Figure 3: Trajectory & Regulon Analysis",
                    theme = theme(plot.title = element_text(size = 14, face = "bold")))
  
  ggsave(file.path(result_dir, "figure3_trajectory_regulon.pdf"), fig3, 
         width = 15, height = 5 * ceiling(length(fig3_plots) / 3))
  ggsave(file.path(result_dir, "figure3_trajectory_regulon.png"), fig3, 
         width = 15, height = 5 * ceiling(length(fig3_plots) / 3), dpi = 200)
}

# =========================================================================
# FIGURE 4: DE & Pathway Summary
# =========================================================================
cat("Generating Figure 4: DE & Pathway Summary...\n")

# Load DE results
if (file.exists("results/05_de/de_tumor_vs_normal_myeloid.csv")) {
  de_tn <- read.csv("results/05_de/de_tumor_vs_normal_myeloid.csv", row.names = 1)
  
  # Mini volcano
  de_tn$significance <- case_when(
    de_tn$p_val_adj < 0.05 & abs(de_tn$avg_log2FC) > 0.5 ~ "Significant",
    TRUE ~ "NS"
  )
  
  # Label top genes
  top_genes <- de_tn %>%
    filter(significance == "Significant") %>%
    arrange(desc(abs(avg_log2FC))) %>%
    head(15)
  
  p4a <- ggplot(de_tn, aes(x = avg_log2FC, y = -log10(p_val_adj), color = significance)) +
    geom_point(size = 0.5, alpha = 0.5) +
    scale_color_manual(values = c("Significant" = "#E74C3C", "NS" = "grey70")) +
    geom_text(data = top_genes, aes(label = gene), size = 2.5,
              nudge_y = 0.5, check_overlap = TRUE) +
    labs(title = "A. Myeloid DE: Tumor vs Normal",
         x = "log2 Fold Change", y = "-log10(adj. p-value)") +
    theme_cowplot() + NoLegend() +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.5) +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", alpha = 0.5)
  
  ggsave(file.path(result_dir, "figure4_de_summary.pdf"), p4a, width = 8, height = 6)
  ggsave(file.path(result_dir, "figure4_de_summary.png"), p4a, width = 8, height = 6, dpi = 200)
}

cat("\n Phase 10 complete! All summary figures saved.\n")
cat(sprintf("Results directory: %s/\n", result_dir))
