# =============================================================================
# 07_trajectory_analysis.R
# Trajectory analysis: Slingshot + Monocle3 on myeloid cells
# Focus: Monocyte-to-macrophage differentiation trajectory
# KEY ADDITION: Compare trajectories across Normal → MMRp → MMRd
# =============================================================================

library(Seurat)
library(slingshot)
library(tradeSeq)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(cowplot)
library(viridis)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

set.seed(42)
cat("=== Phase 7: Trajectory Analysis ===\n\n")

result_dir <- "results/07_trajectory"
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

# Condition colors used throughout
condition_colors <- c("N" = "#3498DB", "MMRp" = "#2ECC71", "MMRd" = "#E74C3C")

# --- 1. Load myeloid Seurat object ---
cat("Loading myeloid Seurat object...\n")
seu <- readRDS("results/04_myeloid/seu_myeloid.rds")

# --- Create a unified condition column: Normal / MMRp / MMRd ---
seu$Condition <- case_when(
  seu$SPECIMEN_TYPE == "N" ~ "N",
  seu$MMRStatus == "MMRp"  ~ "MMRp",
  seu$MMRStatus == "MMRd"  ~ "MMRd",
  TRUE ~ NA_character_
)

# Filter to cells with defined conditions
seu_cond <- subset(seu, subset = !is.na(Condition))
cat(sprintf("  Cells with defined condition (N/MMRp/MMRd): %d\n", ncol(seu_cond)))
cat("  Breakdown:\n")
print(table(seu_cond$Condition))

# =========================================================================
# PART A: SLINGSHOT TRAJECTORY (ALL MYELOID)
# =========================================================================
cat("\n=== Slingshot Trajectory Analysis ===\n")

# --- 2. Run Slingshot on all myeloid cells ---
cat("Preparing data for Slingshot...\n")
umap_coords <- Embeddings(seu, "umap")
clusters <- seu$cl295v11SubShort

cat("Running Slingshot (root: cM01 Monocyte)...\n")
sds <- slingshot(umap_coords, clusterLabels = clusters, start.clus = "cM01")

pseudotime_df <- slingPseudotime(sds)
n_lineages <- ncol(pseudotime_df)
cat(sprintf("  Found %d lineages\n", n_lineages))

for (i in seq_len(n_lineages)) {
  seu@meta.data[[paste0("slingshot_pt", i)]] <- pseudotime_df[, i]
  seu_cond@meta.data[[paste0("slingshot_pt", i)]] <- pseudotime_df[colnames(seu_cond), i]
}

curves <- slingCurves(sds)

# --- 3. Standard Slingshot plots ---
cat("Generating Slingshot trajectory plots...\n")

plot_list <- list()
for (i in seq_len(min(n_lineages, 4))) {
  pt_col <- paste0("slingshot_pt", i)
  df <- data.frame(UMAP_1 = umap_coords[, 1], UMAP_2 = umap_coords[, 2],
                   pseudotime = seu@meta.data[[pt_col]], cluster = clusters)
  curve_df <- as.data.frame(curves[[i]]$s[curves[[i]]$ord, ])
  colnames(curve_df) <- c("UMAP_1", "UMAP_2")
  
  p <- ggplot(df, aes(x = UMAP_1, y = UMAP_2, color = pseudotime)) +
    geom_point(size = 0.3, alpha = 0.5) +
    scale_color_viridis(option = "inferno", na.value = "grey80") +
    geom_path(data = curve_df, aes(x = UMAP_1, y = UMAP_2),
              color = "black", linewidth = 1.5, inherit.aes = FALSE) +
    labs(title = paste0("Slingshot Lineage ", i), color = "Pseudotime") +
    theme_cowplot()
  plot_list[[i]] <- p
}

p_sling <- wrap_plots(plot_list, ncol = 2)
ggsave(file.path(result_dir, "slingshot_pseudotime_umap.pdf"), p_sling,
       width = 12, height = 6 * ceiling(n_lineages / 2))
ggsave(file.path(result_dir, "slingshot_pseudotime_umap.png"), p_sling,
       width = 12, height = 6 * ceiling(n_lineages / 2), dpi = 200)

# Slingshot with cluster labels + all curves
df_clust <- data.frame(UMAP_1 = umap_coords[, 1], UMAP_2 = umap_coords[, 2],
                       cluster = clusters)
p_sling_clust <- ggplot(df_clust, aes(x = UMAP_1, y = UMAP_2, color = cluster)) +
  geom_point(size = 0.3, alpha = 0.5) +
  labs(title = "Slingshot Trajectories on Myeloid UMAP", color = "Cluster") +
  theme_cowplot() +
  guides(color = guide_legend(override.aes = list(size = 3)))
for (i in seq_len(n_lineages)) {
  curve_df <- as.data.frame(curves[[i]]$s[curves[[i]]$ord, ])
  colnames(curve_df) <- c("UMAP_1", "UMAP_2")
  p_sling_clust <- p_sling_clust +
    geom_path(data = curve_df, aes(x = UMAP_1, y = UMAP_2),
              color = "black", linewidth = 1.2, inherit.aes = FALSE)
}
ggsave(file.path(result_dir, "slingshot_trajectories_clusters.pdf"), p_sling_clust, width = 10, height = 7)
ggsave(file.path(result_dir, "slingshot_trajectories_clusters.png"), p_sling_clust, width = 10, height = 7, dpi = 200)

# =========================================================================
# PART B: CONDITION-SPECIFIC TRAJECTORY ANALYSIS
# Two parallel comparisons: Normal → MMRp  and  Normal → MMRd
# =========================================================================
cat("\n=== Condition-Specific Trajectory Analysis ===\n")
cat("  Comparing: (1) Normal → pMMR   (2) Normal → dMMR\n")

# --- 4. UMAP colored by condition with trajectory overlay ---
cat("Plotting trajectory by condition...\n")

umap_cond <- data.frame(
  UMAP_1 = Embeddings(seu_cond, "umap")[, 1],
  UMAP_2 = Embeddings(seu_cond, "umap")[, 2],
  Condition = seu_cond$Condition,
  pseudotime = seu_cond$slingshot_pt1,
  cluster = seu_cond$cl295v11SubShort
)

# UMAP split by condition
# Set condition factor levels (Normal first, then the two tumor types)
umap_cond$Condition <- factor(umap_cond$Condition, levels = c("N", "MMRp", "MMRd"))

p_cond_split <- ggplot(umap_cond, aes(x = UMAP_1, y = UMAP_2, color = pseudotime)) +
  geom_point(size = 0.3, alpha = 0.5) +
  scale_color_viridis(option = "inferno", na.value = "grey80") +
  facet_wrap(~Condition) +
  labs(title = "Pseudotime: Normal vs pMMR (MMRp) vs dMMR (MMRd)") +
  theme_cowplot() +
  theme(strip.text = element_text(size = 12, face = "bold"))

ggsave(file.path(result_dir, "pseudotime_by_condition_umap.pdf"), p_cond_split, width = 16, height = 5)
ggsave(file.path(result_dir, "pseudotime_by_condition_umap.png"), p_cond_split, width = 16, height = 5, dpi = 200)

# UMAP colored by condition (single panel)
p_cond_overlay <- ggplot(umap_cond, aes(x = UMAP_1, y = UMAP_2, color = Condition)) +
  geom_point(size = 0.3, alpha = 0.3) +
  scale_color_manual(values = condition_colors) +
  labs(title = "Myeloid UMAP by Condition") +
  theme_cowplot() +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))

ggsave(file.path(result_dir, "umap_by_condition.pdf"), p_cond_overlay, width = 8, height = 6)
ggsave(file.path(result_dir, "umap_by_condition.png"), p_cond_overlay, width = 8, height = 6, dpi = 200)

# --- 5. SEPARATE pseudotime analyses: Normal→pMMR and Normal→dMMR ---
cat("\n--- Two parallel comparisons: Normal→pMMR and Normal→dMMR ---\n")

# Run separate Slingshot trajectories for each N+tumor pair
for (comparison in list(
  list(name = "N_vs_MMRp", label = "Normal → pMMR", groups = c("N", "MMRp")),
  list(name = "N_vs_MMRd", label = "Normal → dMMR", groups = c("N", "MMRd"))
)) {
  cat(sprintf("\n  === %s ===\n", comparison$label))
  
  # Subset cells
  seu_pair <- subset(seu_cond, subset = Condition %in% comparison$groups)
  cat(sprintf("    Cells: %d (N=%d, %s=%d)\n", ncol(seu_pair),
              sum(seu_pair$Condition == "N"),
              comparison$groups[2], sum(seu_pair$Condition == comparison$groups[2])))
  
  # Run Slingshot on this pair
  pair_umap <- Embeddings(seu_pair, "umap")
  pair_clusters <- seu_pair$cl295v11SubShort
  sds_pair <- slingshot(pair_umap, clusterLabels = pair_clusters, start.clus = "cM01")
  pseudo_pair <- slingPseudotime(sds_pair)
  seu_pair@meta.data$pseudotime_pair <- pseudo_pair[, 1]
  curves_pair <- slingCurves(sds_pair)
  
  # --- Pseudotime density: Normal vs Tumor ---
  pt_dens_df <- data.frame(
    pseudotime = seu_pair$pseudotime_pair,
    Condition = seu_pair$Condition
  ) %>% filter(!is.na(pseudotime))
  
  pair_colors <- c("N" = "#3498DB")
  pair_colors[comparison$groups[2]] <- condition_colors[comparison$groups[2]]
  
  p_dens <- ggplot(pt_dens_df, aes(x = pseudotime, fill = Condition, color = Condition)) +
    geom_density(alpha = 0.3, linewidth = 0.8) +
    scale_fill_manual(values = pair_colors) +
    scale_color_manual(values = pair_colors) +
    labs(title = paste0("Pseudotime Distribution: ", comparison$label),
         x = "Pseudotime", y = "Density") +
    theme_cowplot()
  
  ggsave(file.path(result_dir, paste0("pseudotime_density_", comparison$name, ".pdf")), p_dens, width = 8, height = 5)
  ggsave(file.path(result_dir, paste0("pseudotime_density_", comparison$name, ".png")), p_dens, width = 8, height = 5, dpi = 200)
  
  # --- Violin / box plot ---
  p_box <- ggplot(pt_dens_df, aes(x = Condition, y = pseudotime, fill = Condition)) +
    geom_violin(alpha = 0.6) +
    geom_boxplot(width = 0.15, outlier.size = 0.3, alpha = 0.8) +
    scale_fill_manual(values = pair_colors) +
    labs(title = paste0("Pseudotime: ", comparison$label), y = "Pseudotime") +
    theme_cowplot() + NoLegend()
  
  # Add Wilcoxon p-value
  wt <- wilcox.test(pseudotime ~ Condition, data = pt_dens_df)
  p_box <- p_box + annotate("text", x = 1.5, y = max(pt_dens_df$pseudotime, na.rm = TRUE) * 0.95,
                             label = paste0("Wilcoxon p=", format.pval(wt$p.value, digits = 3)),
                             size = 4, fontface = "italic")
  
  ggsave(file.path(result_dir, paste0("pseudotime_boxplot_", comparison$name, ".pdf")), p_box, width = 6, height = 5)
  ggsave(file.path(result_dir, paste0("pseudotime_boxplot_", comparison$name, ".png")), p_box, width = 6, height = 5, dpi = 200)
  
  # --- UMAP with trajectory colored by condition ---
  pair_df <- data.frame(
    UMAP_1 = pair_umap[, 1], UMAP_2 = pair_umap[, 2],
    Condition = seu_pair$Condition, pseudotime = seu_pair$pseudotime_pair
  )
  curve1_df <- as.data.frame(curves_pair[[1]]$s[curves_pair[[1]]$ord, ])
  colnames(curve1_df) <- c("UMAP_1", "UMAP_2")
  
  p_traj_cond <- ggplot(pair_df, aes(x = UMAP_1, y = UMAP_2, color = Condition)) +
    geom_point(size = 0.3, alpha = 0.4) +
    scale_color_manual(values = pair_colors) +
    geom_path(data = curve1_df, aes(x = UMAP_1, y = UMAP_2),
              color = "black", linewidth = 1.2, inherit.aes = FALSE) +
    labs(title = paste0("Trajectory: ", comparison$label)) +
    theme_cowplot() + guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))
  
  p_traj_pt <- ggplot(pair_df, aes(x = UMAP_1, y = UMAP_2, color = pseudotime)) +
    geom_point(size = 0.3, alpha = 0.5) +
    scale_color_viridis(option = "inferno", na.value = "grey80") +
    geom_path(data = curve1_df, aes(x = UMAP_1, y = UMAP_2),
              color = "black", linewidth = 1.2, inherit.aes = FALSE) +
    labs(title = paste0("Pseudotime: ", comparison$label)) +
    theme_cowplot()
  
  p_pair <- p_traj_cond | p_traj_pt
  ggsave(file.path(result_dir, paste0("trajectory_umap_", comparison$name, ".pdf")), p_pair, width = 14, height = 6)
  ggsave(file.path(result_dir, paste0("trajectory_umap_", comparison$name, ".png")), p_pair, width = 14, height = 6, dpi = 200)
  
  # --- Macrophage subtype proportions along pseudotime ---
  pt_clust <- data.frame(
    pseudotime = seu_pair$pseudotime_pair,
    Condition = seu_pair$Condition,
    Cluster = seu_pair$cl295v11SubFull
  ) %>% filter(!is.na(pseudotime))
  
  pt_clust$pt_bin <- cut(pt_clust$pseudotime,
                          breaks = quantile(pt_clust$pseudotime, probs = seq(0, 1, 0.1)),
                          include.lowest = TRUE, labels = paste0("Q", 1:10))
  
  bin_prop <- pt_clust %>%
    group_by(Condition, pt_bin, Cluster) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(Condition, pt_bin) %>%
    mutate(prop = n / sum(n))
  
  p_stream <- ggplot(bin_prop, aes(x = pt_bin, y = prop, fill = Cluster)) +
    geom_bar(stat = "identity", position = "stack") +
    facet_wrap(~Condition, nrow = 1) +
    labs(title = paste0("Subtype Composition Along Pseudotime: ", comparison$label),
         x = "Pseudotime Bin", y = "Proportion", fill = "Subcluster") +
    theme_cowplot() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
          strip.text = element_text(size = 12, face = "bold"),
          legend.text = element_text(size = 7))
  
  ggsave(file.path(result_dir, paste0("subtype_along_pt_", comparison$name, ".pdf")), p_stream, width = 12, height = 6)
  ggsave(file.path(result_dir, paste0("subtype_along_pt_", comparison$name, ".png")), p_stream, width = 12, height = 6, dpi = 200)
  
  # --- Gene expression dynamics: Normal vs Tumor along pseudotime ---
  key_genes <- c("SPP1", "C1QA", "C1QB", "CD163", "MRC1", "VEGFA",
                 "TNF", "IL1B", "CXCL9", "CXCL10", "IDO1",
                 "S100A8", "S100A9", "APOE", "HIF1A", "CD14",
                 "MARCO", "TREM2", "GPNMB", "FN1")
  key_genes <- key_genes[key_genes %in% rownames(seu_pair)]
  
  DefaultAssay(seu_pair) <- "RNA"
  expr_pair <- GetAssayData(seu_pair, slot = "data", assay = "RNA")
  
  smoothed_plots <- list()
  for (gene in key_genes) {
    gene_df <- data.frame(
      pseudotime = seu_pair$pseudotime_pair,
      expression = as.numeric(expr_pair[gene, ]),
      Condition = seu_pair$Condition
    ) %>% filter(!is.na(pseudotime))
    
    p <- ggplot(gene_df, aes(x = pseudotime, y = expression, color = Condition)) +
      geom_smooth(method = "loess", se = TRUE, alpha = 0.2, span = 0.5, linewidth = 1) +
      scale_color_manual(values = pair_colors) +
      labs(title = gene, x = "Pseudotime", y = "Expression") +
      theme_cowplot() +
      theme(plot.title = element_text(size = 11, face = "bold"),
            legend.position = "none")
    smoothed_plots[[gene]] <- p
  }
  
  p_dynamics <- wrap_plots(smoothed_plots, ncol = 5) +
    plot_annotation(title = paste0("Gene Dynamics Along Pseudotime: ", comparison$label),
                    theme = theme(plot.title = element_text(size = 14, face = "bold")))
  
  ggsave(file.path(result_dir, paste0("gene_dynamics_", comparison$name, ".pdf")), p_dynamics,
         width = 20, height = 4 * ceiling(length(key_genes) / 5))
  ggsave(file.path(result_dir, paste0("gene_dynamics_", comparison$name, ".png")), p_dynamics,
         width = 20, height = 4 * ceiling(length(key_genes) / 5), dpi = 200)
  
  # --- Condition-specific trajectory heatmaps ---
  for (cond in comparison$groups) {
    cond_cells <- colnames(seu_pair)[seu_pair$Condition == cond]
    cond_pt <- seu_pair$pseudotime_pair[cond_cells]
    valid <- !is.na(cond_pt)
    cond_cells <- cond_cells[valid]
    cond_pt <- cond_pt[valid]
    if (length(cond_cells) < 50) next
    
    n_bins <- 50
    bin_breaks <- seq(min(cond_pt), max(cond_pt), length.out = n_bins + 1)
    bin_lab <- cut(cond_pt, breaks = bin_breaks, include.lowest = TRUE)
    cond_expr <- expr_pair[key_genes, cond_cells]
    
    binned <- sapply(levels(bin_lab), function(b) {
      cells <- which(bin_lab == b)
      if (length(cells) > 0) Matrix::rowMeans(cond_expr[, cells, drop = FALSE])
      else rep(NA, length(key_genes))
    })
    rownames(binned) <- key_genes
    binned <- binned[, !apply(binned, 2, function(x) all(is.na(x)))]
    binned_scaled <- t(scale(t(binned)))
    binned_scaled[binned_scaled > 2] <- 2
    binned_scaled[binned_scaled < -2] <- -2
    
    cond_label <- ifelse(cond == "N", "Normal", cond)
    pdf(file.path(result_dir, paste0("trajectory_heatmap_", comparison$name, "_", cond, ".pdf")), width = 10, height = 6)
    ht <- Heatmap(binned_scaled, name = "Z-score",
                  col = colorRamp2(c(-2, 0, 2), c("#3498DB", "white", "#E74C3C")),
                  cluster_columns = FALSE, cluster_rows = TRUE,
                  show_column_names = FALSE,
                  column_title = paste0(cond_label, " (", comparison$label, "): Gene Expression Along Pseudotime →"),
                  row_names_gp = gpar(fontsize = 8))
    draw(ht)
    dev.off()
  }
  
  cat(sprintf("  Completed: %s\n", comparison$label))
}

# --- 8. Heatmap: condition-specific gene expression along pseudotime ---
cat("Generating condition-specific trajectory heatmaps...\n")

expr_data <- GetAssayData(seu_cond, slot = "data", assay = "RNA")

for (cond in c("N", "MMRp", "MMRd")) {
  cond_cells <- colnames(seu_cond)[seu_cond$Condition == cond]
  cond_pt <- seu_cond$slingshot_pt1[cond_cells]
  valid <- !is.na(cond_pt)
  cond_cells <- cond_cells[valid]
  cond_pt <- cond_pt[valid]
  
  if (length(cond_cells) < 50) next
  
  # Bin into 50 bins
  n_bins <- 50
  bin_breaks <- seq(min(cond_pt), max(cond_pt), length.out = n_bins + 1)
  bin_labels <- cut(cond_pt, breaks = bin_breaks, include.lowest = TRUE)
  
  cond_expr <- expr_data[key_genes, cond_cells]
  
  binned <- sapply(levels(bin_labels), function(b) {
    cells <- which(bin_labels == b)
    if (length(cells) > 0) Matrix::rowMeans(cond_expr[, cells, drop = FALSE])
    else rep(NA, length(key_genes))
  })
  rownames(binned) <- key_genes
  binned <- binned[, !apply(binned, 2, function(x) all(is.na(x)))]
  binned_scaled <- t(scale(t(binned)))
  binned_scaled[binned_scaled > 2] <- 2
  binned_scaled[binned_scaled < -2] <- -2
  
  pdf(file.path(result_dir, paste0("trajectory_heatmap_", cond, ".pdf")), width = 10, height = 6)
  ht <- Heatmap(binned_scaled, name = "Z-score",
                col = colorRamp2(c(-2, 0, 2), c("#3498DB", "white", "#E74C3C")),
                cluster_columns = FALSE, cluster_rows = TRUE,
                show_column_names = FALSE,
                column_title = paste0(cond, ": Gene Expression Along Pseudotime →"),
                row_names_gp = gpar(fontsize = 8))
  draw(ht)
  dev.off()
}

# =========================================================================
# PART C: tradeSeq - trajectory-associated genes
# =========================================================================
cat("\n=== tradeSeq: Trajectory-Associated Genes ===\n")

counts_mat <- GetAssayData(seu, slot = "counts", assay = "RNA")
n_cells <- ncol(counts_mat)

if (n_cells > 10000) {
  cat(sprintf("  Downsampling from %d to 10000 cells for tradeSeq...\n", n_cells))
  set.seed(42)
  sample_idx <- sample(n_cells, 10000)
  counts_sub <- counts_mat[, sample_idx]
  pseudo_sub <- pseudotime_df[sample_idx, , drop = FALSE]
  weights_sub <- slingCurveWeights(sds)[sample_idx, , drop = FALSE]
} else {
  counts_sub <- counts_mat
  pseudo_sub <- pseudotime_df
  weights_sub <- slingCurveWeights(sds)
}

gene_filter <- rowSums(counts_sub > 0) >= 50
counts_sub <- counts_sub[gene_filter, ]
cat(sprintf("  Using %d genes, %d cells\n", nrow(counts_sub), ncol(counts_sub)))

cat("  Fitting GAM...\n")
sce <- fitGAM(counts = as.matrix(counts_sub), pseudotime = pseudo_sub,
              cellWeights = weights_sub, nknots = 6, verbose = TRUE)

cat("  Running association test...\n")
assoc_res <- associationTest(sce, lineages = TRUE)
assoc_res$gene <- rownames(assoc_res)
assoc_res <- assoc_res %>% arrange(pvalue)

write.csv(assoc_res, file.path(result_dir, "tradeseq_association_results.csv"), row.names = FALSE)

top_traj_genes <- head(assoc_res$gene[assoc_res$pvalue < 0.01], 50)
cat(sprintf("  Trajectory-associated genes (p < 0.01): %d\n",
            sum(assoc_res$pvalue < 0.01, na.rm = TRUE)))

# Top trajectory genes heatmap (all cells)
if (length(top_traj_genes) > 5) {
  cat("Generating trajectory gene heatmap (all cells)...\n")
  
  heatmap_genes <- intersect(top_traj_genes, rownames(seu))
  valid_cells <- !is.na(seu$slingshot_pt1)
  valid_pt <- seu$slingshot_pt1[valid_cells]
  valid_expr <- GetAssayData(seu, slot = "data", assay = "RNA")[heatmap_genes, valid_cells]
  
  n_bins <- 100
  bin_breaks <- seq(min(valid_pt), max(valid_pt), length.out = n_bins + 1)
  bin_labels <- cut(valid_pt, breaks = bin_breaks, include.lowest = TRUE)
  
  binned_expr <- sapply(levels(bin_labels), function(b) {
    cells <- which(bin_labels == b)
    if (length(cells) > 0) Matrix::rowMeans(valid_expr[, cells, drop = FALSE])
    else rep(NA, length(heatmap_genes))
  })
  rownames(binned_expr) <- heatmap_genes
  binned_expr <- binned_expr[, !apply(binned_expr, 2, function(x) all(is.na(x)))]
  binned_scaled <- t(scale(t(binned_expr)))
  binned_scaled[binned_scaled > 2] <- 2
  binned_scaled[binned_scaled < -2] <- -2
  
  pdf(file.path(result_dir, "trajectory_gene_heatmap.pdf"), width = 12, height = 10)
  ht <- Heatmap(binned_scaled, name = "Z-score",
                col = colorRamp2(c(-2, 0, 2), c("#3498DB", "white", "#E74C3C")),
                cluster_columns = FALSE, cluster_rows = TRUE,
                show_column_names = FALSE,
                column_title = "Pseudotime (Lineage 1) →",
                row_names_gp = gpar(fontsize = 7), use_raster = TRUE)
  draw(ht)
  dev.off()
}

# =========================================================================
# PART D: MONOCLE3 (if available)
# =========================================================================
cat("\n=== Monocle3 Trajectory Analysis ===\n")

if (requireNamespace("monocle3", quietly = TRUE) &
    requireNamespace("SeuratWrappers", quietly = TRUE)) {
  library(monocle3)
  library(SeuratWrappers)
  
  cat("Converting Seurat to Monocle3 CDS...\n")
  cds <- as.cell_data_set(seu, assay = "RNA")
  # Inject existing Seurat UMAP coords and Harmony-corrected cluster labels.
  # Avoids re-running Monocle3's own clustering, which would produce different
  # partitions and break cross-method comparability with Slingshot (Part A/B).
  recreate_partitions <- rep(1, ncol(cds))
  names(recreate_partitions) <- colnames(cds)
  cds@clusters@listData[["UMAP"]][["partitions"]] <- as.factor(recreate_partitions)
  cds@clusters@listData[["UMAP"]][["clusters"]]   <- seu$cl295v11SubShort[colnames(cds)]
  cds@int_colData@listData[["reducedDims"]]@listData[["UMAP"]] <-
    Embeddings(seu, "umap")[colnames(cds), ]

  cat("Learning trajectory graph...\n")
  cds <- learn_graph(cds, use_partition = FALSE)

  cat("Ordering cells (root: Monocytes)...\n")
  monocyte_cells <- colnames(cds)[seu$cl295v11SubShort[colnames(cds)] == "cM01"]
  cds <- order_cells(cds, root_cells = monocyte_cells[1:min(10, length(monocyte_cells))])
  
  seu$monocle3_pseudotime <- pseudotime(cds)
  
  # Standard plots
  p_m3_pt <- plot_cells(cds, color_cells_by = "pseudotime",
                         label_cell_groups = FALSE, label_leaves = FALSE,
                         label_branch_points = FALSE, cell_size = 0.3) +
    scale_color_viridis(option = "inferno") +
    labs(title = "Monocle3 Pseudotime") + theme_cowplot()
  
  p_m3_cl <- plot_cells(cds, color_cells_by = "cl295v11SubFull",
                         label_cell_groups = TRUE, label_leaves = FALSE,
                         label_branch_points = FALSE, cell_size = 0.3) +
    labs(title = "Monocle3 Trajectory - Clusters") + theme_cowplot()
  
  # Condition-specific Monocle3 UMAP
  # Add condition to CDS
  colData(cds)$Condition <- seu$Condition[colnames(cds)]
  
  p_m3_cond <- plot_cells(cds, color_cells_by = "Condition",
                           label_cell_groups = FALSE, label_leaves = FALSE,
                           label_branch_points = FALSE, cell_size = 0.3) +
    scale_color_manual(values = condition_colors, na.value = "grey80") +
    labs(title = "Monocle3 - By Condition") + theme_cowplot()
  
  p_m3_all <- (p_m3_cl | p_m3_pt | p_m3_cond)
  ggsave(file.path(result_dir, "monocle3_trajectory.pdf"), p_m3_all, width = 21, height = 7)
  ggsave(file.path(result_dir, "monocle3_trajectory.png"), p_m3_all, width = 21, height = 7, dpi = 200)
  
  # Monocle3 pseudotime density by condition
  m3_pt_data <- data.frame(
    pseudotime = pseudotime(cds),
    Condition = colData(cds)$Condition
  ) %>% filter(!is.na(Condition) & is.finite(pseudotime))
  
  p_m3_dens <- ggplot(m3_pt_data, aes(x = pseudotime, fill = Condition, color = Condition)) +
    geom_density(alpha = 0.3, linewidth = 0.8) +
    scale_fill_manual(values = condition_colors) +
    scale_color_manual(values = condition_colors) +
    labs(title = "Monocle3 Pseudotime Distribution: N vs MMRp vs MMRd",
         x = "Pseudotime", y = "Density") +
    theme_cowplot()
  
  ggsave(file.path(result_dir, "monocle3_pseudotime_density_by_condition.pdf"), p_m3_dens, width = 8, height = 5)
  ggsave(file.path(result_dir, "monocle3_pseudotime_density_by_condition.png"), p_m3_dens, width = 8, height = 5, dpi = 200)
  
  # Graph test
  cat("Running Monocle3 graph_test...\n")
  graph_test_res <- graph_test(cds, neighbor_graph = "principal_graph", cores = 4)
  graph_test_res <- graph_test_res %>% arrange(q_value)
  write.csv(graph_test_res, file.path(result_dir, "monocle3_graph_test_results.csv"), row.names = FALSE)
  cat(sprintf("  Trajectory-varying genes (q < 0.05): %d\n",
              sum(graph_test_res$q_value < 0.05, na.rm = TRUE)))
  
} else {
  cat("  monocle3 or SeuratWrappers not available. Skipping.\n")
  cat("  Install: devtools::install_github('cole-trapnell-lab/monocle3')\n")
  cat("  Also: remotes::install_github('satijalab/seurat-wrappers')\n")
}

# =========================================================================
# PART E: DIFFUSIONMAP-BASED CELL PROPENSITY SCORES
# Adapted from Karimi et al. (eLife 2024, MAM macrophage paper):
# Quantifies P(transitional myeloid cell reaches a terminal macrophage state)
# per condition, revealing condition-specific fate commitment beyond pseudotime.
# Dependencies: destiny, ggpubr
# =========================================================================
cat("\n=== DiffusionMap Propensity Scoring ===\n")

if (requireNamespace("destiny",  quietly = TRUE) &&
    requireNamespace("ggpubr",   quietly = TRUE)) {
  library(destiny)
  library(ggpubr)

  # TODO(human): Define cluster roles for the CRC myeloid trajectory.
  # Look at the marker outputs from 04_myeloid_analysis.R (dot plots, FindAllMarkers)
  # to classify each cM01–cM10 cluster as:
  #   "source"      – monocytes, trajectory root (cM01)
  #   "transitional" – intermediate clusters without a committed phenotype
  #   "terminal"    – terminally polarised macrophage end-states
  # Also identify which terminal clusters are enriched in MMRp vs MMRd tumors
  # (check the composition plots from 03_composition.R or 04_myeloid_analysis.R).
  cluster_roles <- list(
    transitional  = c("cM01"),       # TODO(human): add intermediate clusters
    terminal_mmrp = c("cM02"),       # TODO(human): MMRp-enriched terminal clusters
    terminal_mmrd = c("cM02")        # TODO(human): MMRd-enriched terminal clusters
  )

  propensity_calc <- function(seu_sub, dm, trans_clust, recv_clust, targ_clust,
                               p_thresh = 0.05) {
    Trans_t  <- dm@transitions / rowSums(dm@transitions)
    int_c    <- intersect(colnames(Trans_t), colnames(seu_sub))
    Trans_t  <- Trans_t[int_c, int_c]
    is_valid <- rowSums(is.na(Trans_t)) == 0
    clust    <- seu_sub$cl295v11SubShort[int_c]
    cond     <- seu_sub$Condition[int_c]

    do.call(rbind, lapply(unique(na.omit(cond)), function(d) {
      is_d   <- cond == d & !is.na(cond)
      i_src  <- clust %in% trans_clust & is_valid & is_d
      i_recv <- clust %in% recv_clust  & is_valid & is_d
      i_targ <- clust %in% targ_clust  & is_valid & is_d
      if (sum(i_src) < 2) return(NULL)

      r <- if (sum(i_recv) > 1) rowSums(Trans_t[i_src, i_recv, drop = FALSE]) else
               Trans_t[i_src, which(i_recv), drop = TRUE]
      s <- if (sum(i_targ) > 1) rowSums(Trans_t[i_src, i_targ, drop = FALSE]) else
               Trans_t[i_src, which(i_targ), drop = TRUE]

      q <- r + s; keep <- q > p_thresh
      if (!any(keep)) return(NULL)
      data.frame(Condition = d, propensity = (s / (r + s))[keep],
                 cell_id = rownames(Trans_t)[i_src][keep])
    }))
  }

  all_relevant <- unique(c(cluster_roles$transitional,
                            cluster_roles$terminal_mmrp,
                            cluster_roles$terminal_mmrd))
  seu_dm <- subset(seu_cond, subset = cl295v11SubShort %in% all_relevant)
  cat(sprintf("  DiffusionMap input: %d cells\n", ncol(seu_dm)))

  dm_path <- file.path(result_dir, "diffusionmap_myeloid.rds")
  if (file.exists(dm_path)) {
    cat("  Loading saved DiffusionMap...\n")
    dm <- readRDS(dm_path)
  } else {
    cat("  Building DiffusionMap (saved to rds for reuse)...\n")
    sce_dm <- as.SingleCellExperiment(seu_dm, assay = "RNA")
    dm <- DiffusionMap(sce_dm, verbose = FALSE,
                       k = round(0.05 * ncol(seu_dm)), n_pcs = 30)
    saveRDS(dm, dm_path)
  }

  prop_df <- propensity_calc(
    seu_sub    = seu_dm,
    dm         = dm,
    trans_clust = cluster_roles$transitional,
    recv_clust  = cluster_roles$terminal_mmrp,
    targ_clust  = cluster_roles$terminal_mmrd
  )

  if (!is.null(prop_df) && nrow(prop_df) > 0) {
    prop_df$Condition <- factor(prop_df$Condition, levels = c("N", "MMRp", "MMRd"))

    p_prop_vln <- ggplot(prop_df,
                         aes(x = Condition, y = propensity, fill = Condition)) +
      geom_violin(adjust = 2, scale = "width", alpha = 0.8) +
      geom_jitter(height = 0, width = 0.15, size = 0.5, alpha = 0.4) +
      stat_summary(fun = median, geom = "crossbar",
                   width = 0.5, colour = "white", linewidth = 0.3) +
      scale_fill_manual(values = condition_colors) +
      stat_compare_means(comparisons = list(c("MMRp", "N"), c("MMRd", "N")),
                         method = "wilcox.test", size = 3.5) +
      labs(title = "Macrophage Differentiation Propensity by Condition",
           subtitle = "P(transitional cell reaches terminal macrophage state)",
           y = "Propensity score", x = "") +
      theme_cowplot() + NoLegend()

    ggsave(file.path(result_dir, "propensity_violin_by_condition.pdf"),
           p_prop_vln, width = 6, height = 6)
    ggsave(file.path(result_dir, "propensity_violin_by_condition.png"),
           p_prop_vln, width = 6, height = 6, dpi = 200)
    write.csv(prop_df, file.path(result_dir, "propensity_scores_by_condition.csv"),
              row.names = FALSE)
    cat("  Propensity scores saved.\n")
  }

} else {
  cat("  destiny or ggpubr not available. Skipping propensity analysis.\n")
  cat("  Install: install.packages(c('destiny', 'ggpubr'))\n")
}

# --- Save ---
cat("\nSaving myeloid object with pseudotime + condition...\n")
saveRDS(seu, file.path(result_dir, "seu_myeloid_trajectory.rds"))

cat("\n Phase 7 complete!\n")

# Helper function (defined at end to avoid errors if package not available)
stat_compare_means_if_available <- function(data) {
  # Simple annotation of Kruskal-Wallis p-value
  kw <- kruskal.test(pseudotime ~ Condition, data = data)
  annotate("text", x = 2, y = max(data$pseudotime, na.rm = TRUE) * 0.95,
           label = paste0("Kruskal-Wallis p=", format.pval(kw$p.value, digits = 3)),
           size = 3.5, fontface = "italic")
}
