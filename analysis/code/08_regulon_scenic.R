# =============================================================================
# 08_regulon_scenic.R
# Regulon analysis using SCENIC on myeloid cells
# KEY ADDITION: Compare regulon activity across Normal → MMRp → MMRd
# =============================================================================

library(Seurat)
library(SCENIC)
library(AUCell)
library(RcisTarget)
library(GENIE3)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)
library(cowplot)
library(viridis)
library(patchwork)

set.seed(42)
cat("=== Phase 8: Regulon Analysis (SCENIC) ===\n\n")

result_dir <- "results/08_scenic"
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

# Condition colors
condition_colors <- c("N" = "#3498DB", "MMRp" = "#2ECC71", "MMRd" = "#E74C3C")

# --- 1. Load myeloid Seurat object ---
cat("Loading myeloid Seurat object...\n")
seu <- readRDS("results/04_myeloid/seu_myeloid.rds")

# Create condition column
seu$Condition <- case_when(
  seu$SPECIMEN_TYPE == "N" ~ "N",
  seu$MMRStatus == "MMRp"  ~ "MMRp",
  seu$MMRStatus == "MMRd"  ~ "MMRd",
  TRUE ~ NA_character_
)

cat(sprintf("  Cells with condition: N=%d, MMRp=%d, MMRd=%d\n",
            sum(seu$Condition == "N", na.rm = TRUE),
            sum(seu$Condition == "MMRp", na.rm = TRUE),
            sum(seu$Condition == "MMRd", na.rm = TRUE)))

# =========================================================================
# PART A: Prepare expression matrix
# =========================================================================
cat("\nPreparing expression matrix for SCENIC...\n")
DefaultAssay(seu) <- "RNA"
expr_mat <- GetAssayData(seu, slot = "data", assay = "RNA")

# Downsample - stratified by cluster AND condition
n_cells <- ncol(expr_mat)
if (n_cells > 5000) {
  cat(sprintf("  Downsampling from %d to ~5000 cells (stratified by cluster + condition)...\n", n_cells))
  set.seed(42)
  meta <- seu@meta.data
  meta$cell_id <- rownames(meta)
  
  sampled_cells <- meta %>%
    filter(!is.na(Condition)) %>%
    group_by(cl295v11SubShort, Condition) %>%
    slice_sample(n = min(n(), 200)) %>%
    ungroup() %>%
    pull(cell_id)
  
  if (length(sampled_cells) > 6000) {
    sampled_cells <- sample(sampled_cells, 6000)
  }
  
  expr_mat_sub <- as.matrix(expr_mat[, sampled_cells])
  meta_sub <- seu@meta.data[sampled_cells, ]
} else {
  expr_mat_sub <- as.matrix(expr_mat)
  meta_sub <- seu@meta.data
}

gene_filter <- rowSums(expr_mat_sub > 0) >= ncol(expr_mat_sub) * 0.01
expr_mat_sub <- expr_mat_sub[gene_filter, ]
cat(sprintf("  Matrix: %d genes x %d cells\n", nrow(expr_mat_sub), ncol(expr_mat_sub)))

# =========================================================================
# PART B: cisTarget databases check + SCENIC / curated approach
# =========================================================================
db_dir <- file.path(result_dir, "cisTarget_databases")
dir.create(db_dir, recursive = TRUE, showWarnings = FALSE)

db_files <- c(
  "hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.genes_vs_motifs.rankings.feather",
  "hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather"
)
db_available <- all(file.exists(file.path(db_dir, db_files)))

if (!db_available) {
  cat("\n  cisTarget databases not found. Using curated regulon approach.\n")
  cat("  For full SCENIC, download from: https://resources.aertslab.org/cistarget/\n\n")
}

# =========================================================================
# PART C: Define curated myeloid TF regulons
# =========================================================================
cat("Defining curated myeloid TF regulons...\n")

myeloid_regulons <- list(
  # Core myeloid differentiation
  "SPI1_PU1"  = c("SPI1", "CSF1R", "CD68", "AIF1", "FCGR1A", "FCGR3A", "CD14", "TYROBP", "FCER1G", "LYZ"),
  "CEBPB"     = c("CEBPB", "IL6", "IL1B", "TNF", "S100A8", "S100A9", "IL8", "CCL3", "CCL4", "PTGS2"),
  "CEBPA"     = c("CEBPA", "CSF3R", "MPO", "ELANE", "PRTN3", "AZU1"),
  
  # M2 / tissue-resident / anti-inflammatory
  "MAFB"      = c("MAFB", "APOE", "C1QA", "C1QB", "C1QC", "CD163", "MRC1", "STAB1", "FOLR2", "LYVE1"),
  "MAF"       = c("MAF", "CCL18", "CD163", "MRC1", "IL10", "TGFBI", "MSR1", "SLC40A1", "HMOX1", "SELENOP"),
  "PPARG"     = c("PPARG", "CD36", "FABP4", "LPL", "MSR1", "CD163", "MRC1", "SCARB1", "ABCA1", "HMOX1"),
  "NR1H3_LXR" = c("NR1H3", "ABCA1", "ABCG1", "APOE", "LPL", "SREBF1", "SCD", "FASN", "APOC1"),
  
  # M1 / pro-inflammatory / IFN-response
  "STAT1"     = c("STAT1", "IRF1", "GBP1", "GBP2", "GBP4", "GBP5", "CXCL9", "CXCL10", "CXCL11", "IDO1", "TAP1", "PSMB8"),
  "IRF1"      = c("IRF1", "GBP1", "GBP2", "TAP1", "PSMB9", "STAT1", "CXCL9", "CXCL10", "HLA-A", "HLA-B"),
  "NFkB_RELA" = c("RELA", "NFKBIA", "NFKB1", "TNF", "IL1B", "IL6", "CCL2", "CCL3", "CXCL8", "PTGS2", "SOD2"),
  "STAT3"     = c("STAT3", "BCL2", "MCL1", "MYC", "CCND1", "VEGFA", "HIF1A", "IL10", "SOCS3"),
  
  # Hypoxia / pro-tumoral macrophage
  "HIF1A"     = c("HIF1A", "VEGFA", "SPP1", "LDHA", "SLC2A1", "PGK1", "ENO1", "PKM", "BNIP3", "PDK1", "HILPDA"),
  "MYC"       = c("MYC", "LDHA", "PKM", "NPM1", "NCL", "GAPDH"),
  
  # DC-associated
  "IRF8"      = c("IRF8", "BATF3", "CLEC9A", "XCR1", "IDO1", "THBD", "CADM1", "SNX22", "WDFY4"),
  "IRF4"      = c("IRF4", "CD1C", "FCER1A", "CLEC10A", "CD1E", "HLA-DQA1", "HLA-DQB1"),
  "IRF7"      = c("IRF7", "LILRA4", "IL3RA", "JCHAIN", "GZMB", "TCF4", "ITM2C", "PLD4"),
  
  # Lipid-associated / TREM2+ macrophage
  "TREM2"     = c("TREM2", "GPNMB", "FABP5", "SPP1", "LGALS3", "CD9", "LIPA", "APOE", "CTSD"),
  
  # Tissue remodeling
  "EGR1"      = c("EGR1", "FOS", "JUN", "JUNB", "ATF3", "DUSP1", "ZFP36", "NR4A1"),
  "NFAT"      = c("NFATC1", "NFATC2", "IL2", "TNF", "CSF2", "IFNG")
)

# Filter to genes present
myeloid_regulons_filtered <- lapply(myeloid_regulons, function(genes) {
  genes[genes %in% rownames(expr_mat_sub)]
})
myeloid_regulons_filtered <- myeloid_regulons_filtered[sapply(myeloid_regulons_filtered, length) >= 3]
cat(sprintf("  %d regulons with ≥3 genes\n", length(myeloid_regulons_filtered)))

# =========================================================================
# PART D: Run AUCell
# =========================================================================
cat("\nRunning AUCell...\n")
cells_rankings <- AUCell_buildRankings(expr_mat_sub, plotStats = FALSE, verbose = FALSE)
cells_AUC <- AUCell_calcAUC(myeloid_regulons_filtered, cells_rankings, verbose = FALSE)
regulonAUC <- getAUC(cells_AUC)

cat(sprintf("  AUC matrix: %d regulons x %d cells\n", nrow(regulonAUC), ncol(regulonAUC)))

# =========================================================================
# PART E: Regulon activity heatmap by cluster
# =========================================================================
cat("\nGenerating regulon activity heatmap by cluster...\n")

cluster_labels <- meta_sub$cl295v11SubFull
avg_auc <- sapply(unique(cluster_labels), function(cl) {
  cells <- which(cluster_labels == cl)
  rowMeans(regulonAUC[, cells, drop = FALSE])
})

avg_auc_scaled <- t(scale(t(avg_auc)))
avg_auc_scaled[avg_auc_scaled > 2] <- 2
avg_auc_scaled[avg_auc_scaled < -2] <- -2

pdf(file.path(result_dir, "regulon_activity_heatmap.pdf"), width = 10, height = 10)
ht <- Heatmap(avg_auc_scaled, name = "Scaled\nAUC",
              col = colorRamp2(c(-2, 0, 2), c("#3498DB", "white", "#E74C3C")),
              row_names_gp = gpar(fontsize = 9),
              column_names_gp = gpar(fontsize = 8),
              column_names_rot = 45,
              cluster_columns = TRUE, cluster_rows = TRUE,
              column_title = "Regulon Activity by Myeloid Subcluster",
              rect_gp = gpar(col = "grey80", lwd = 0.5))
draw(ht)
dev.off()

# =========================================================================
# PART F: CONDITION-SPECIFIC REGULON ANALYSIS (N vs MMRp vs MMRd)
# =========================================================================
cat("\n=== Condition-Specific Regulon Analysis ===\n")

# --- F1. Regulon activity heatmap: cluster × condition ---
cat("Computing regulon activity by cluster × condition...\n")

meta_sub$cluster_condition <- paste0(meta_sub$cl295v11SubFull, " | ", meta_sub$Condition)
valid_rows <- !is.na(meta_sub$Condition)

avg_auc_cc <- sapply(unique(meta_sub$cluster_condition[valid_rows]), function(cc) {
  cells <- which(meta_sub$cluster_condition == cc & valid_rows)
  if (length(cells) >= 5) rowMeans(regulonAUC[, cells, drop = FALSE])
  else rep(NA, nrow(regulonAUC))
})

# Remove groups with too few cells
avg_auc_cc <- avg_auc_cc[, !apply(avg_auc_cc, 2, function(x) any(is.na(x)))]
rownames(avg_auc_cc) <- rownames(regulonAUC)

avg_auc_cc_scaled <- t(scale(t(avg_auc_cc)))
avg_auc_cc_scaled[avg_auc_cc_scaled > 2] <- 2
avg_auc_cc_scaled[avg_auc_cc_scaled < -2] <- -2

# Annotation: extract condition from column name
col_conditions <- gsub(".*\\| ", "", colnames(avg_auc_cc_scaled))
col_clusters <- gsub(" \\|.*", "", colnames(avg_auc_cc_scaled))

ha <- HeatmapAnnotation(
  Condition = col_conditions,
  Cluster = col_clusters,
  col = list(
    Condition = condition_colors
  ),
  annotation_name_side = "left"
)

pdf(file.path(result_dir, "regulon_heatmap_cluster_x_condition.pdf"), width = 16, height = 10)
ht <- Heatmap(avg_auc_cc_scaled, name = "Scaled\nAUC",
              col = colorRamp2(c(-2, 0, 2), c("#3498DB", "white", "#E74C3C")),
              top_annotation = ha,
              row_names_gp = gpar(fontsize = 8),
              column_names_gp = gpar(fontsize = 6),
              column_names_rot = 60,
              cluster_columns = TRUE, cluster_rows = TRUE,
              column_title = "Regulon Activity: Cluster × Condition\n(Normal vs pMMR vs dMMR)")
draw(ht)
dev.off()

# --- F2. Regulon activity by condition WITHIN macrophages (cM02) ---
cat("Comparing regulon activity in macrophages (cM02) across conditions...\n")

macro_cells <- which(meta_sub$cl295v11SubShort == "cM02" & !is.na(meta_sub$Condition))

if (length(macro_cells) > 50) {
  avg_auc_macro_cond <- sapply(c("N", "MMRp", "MMRd"), function(cond) {
    cells <- which(meta_sub$cl295v11SubShort == "cM02" & meta_sub$Condition == cond)
    if (length(cells) >= 5) rowMeans(regulonAUC[, cells, drop = FALSE])
    else rep(NA, nrow(regulonAUC))
  })
  rownames(avg_auc_macro_cond) <- rownames(regulonAUC)
  
  # Remove NAs
  valid_regs <- !apply(avg_auc_macro_cond, 1, function(x) any(is.na(x)))
  avg_auc_macro_cond <- avg_auc_macro_cond[valid_regs, ]
  
  # Compute changes relative to Normal for each tumor type separately
  fc_mmrp <- avg_auc_macro_cond[, "MMRp"] - avg_auc_macro_cond[, "N"]
  fc_mmrd <- avg_auc_macro_cond[, "MMRd"] - avg_auc_macro_cond[, "N"]
  
  fc_df <- data.frame(
    Regulon = rownames(avg_auc_macro_cond),
    AUC_Normal = avg_auc_macro_cond[, "N"],
    AUC_pMMR = avg_auc_macro_cond[, "MMRp"],
    AUC_dMMR = avg_auc_macro_cond[, "MMRd"],
    Delta_pMMR_vs_Normal = fc_mmrp,
    Delta_dMMR_vs_Normal = fc_mmrd
  ) %>% arrange(desc(abs(Delta_dMMR_vs_Normal)))
  
  write.csv(fc_df, file.path(result_dir, "regulon_macrophage_condition_comparison.csv"), row.names = FALSE)
  
  cat("\n--- Regulon Activity Changes in Macrophages (vs Normal) ---\n")
  cat("  Two parallel comparisons: Normal→pMMR  and  Normal→dMMR\n")
  print(head(fc_df[, c("Regulon", "Delta_pMMR_vs_Normal", "Delta_dMMR_vs_Normal")], 10))
  
  # Bar plot: regulon changes
  fc_long <- fc_df %>%
    select(Regulon, Delta_pMMR_vs_Normal, Delta_dMMR_vs_Normal) %>%
    pivot_longer(cols = c(Delta_pMMR_vs_Normal, Delta_dMMR_vs_Normal),
                 names_to = "Comparison", values_to = "Delta_AUC") %>%
    mutate(Comparison = case_when(
      Comparison == "Delta_pMMR_vs_Normal" ~ "Normal → pMMR",
      Comparison == "Delta_dMMR_vs_Normal" ~ "Normal → dMMR"
    ))
  
  # Order by MMRd change
  reg_order <- fc_df %>% arrange(Delta_dMMR_vs_Normal) %>% pull(Regulon)
  fc_long$Regulon <- factor(fc_long$Regulon, levels = reg_order)
  
  p_bar <- ggplot(fc_long, aes(x = Regulon, y = Delta_AUC, fill = Comparison)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_manual(values = c("Normal → pMMR" = "#2ECC71", "Normal → dMMR" = "#E74C3C")) +
    coord_flip() +
    labs(title = "Regulon Activity Changes in Macrophages (cM02)",
         subtitle = "Two parallel comparisons from Normal tissue",
         x = "", y = "ΔAUC (vs Normal)") +
    theme_cowplot() +
    geom_hline(yintercept = 0, linetype = "dashed")
  
  ggsave(file.path(result_dir, "regulon_condition_barplot_macrophages.pdf"), p_bar, width = 10, height = 8)
  ggsave(file.path(result_dir, "regulon_condition_barplot_macrophages.png"), p_bar, width = 10, height = 8, dpi = 200)
  
  # Heatmap: macrophage regulon activity by condition
  avg_scaled <- t(scale(t(avg_auc_macro_cond)))
  avg_scaled[avg_scaled > 2] <- 2
  avg_scaled[avg_scaled < -2] <- -2
  
  ha_cond <- HeatmapAnnotation(
    Condition = c("N", "MMRp", "MMRd"),
    col = list(Condition = condition_colors)
  )
  
  pdf(file.path(result_dir, "regulon_heatmap_macrophage_by_condition.pdf"), width = 6, height = 10)
  ht <- Heatmap(avg_scaled, name = "Scaled\nAUC",
                col = colorRamp2(c(-2, 0, 2), c("#3498DB", "white", "#E74C3C")),
                top_annotation = ha_cond,
                row_names_gp = gpar(fontsize = 9),
                column_names_gp = gpar(fontsize = 11),
                cluster_columns = FALSE,
                column_order = c("N", "MMRp", "MMRd"),
                cluster_rows = TRUE,
                column_title = "Macrophage (cM02) Regulon Activity\nNormal vs pMMR vs dMMR (parallel comparisons)")
  draw(ht)
  dev.off()
}

# --- F3. Violin plots of key regulons by condition (all myeloid) ---
cat("\nGenerating regulon violin plots by condition...\n")

# Score regulons on full dataset
for (reg_name in names(myeloid_regulons_filtered)) {
  clean_name <- gsub("[^[:alnum:]]", "_", reg_name)
  seu <- AddModuleScore(seu, features = list(myeloid_regulons_filtered[[reg_name]]),
                        name = paste0("Reg_", clean_name))
}

# Focus on key regulons
key_regulons <- c("STAT1", "HIF1A", "MAFB", "CEBPB", "NFkB_RELA", "TREM2",
                  "SPI1_PU1", "IRF1", "MAF", "PPARG")
key_regulons <- key_regulons[key_regulons %in% names(myeloid_regulons_filtered)]

violin_plots <- list()
for (reg in key_regulons) {
  clean_name <- gsub("[^[:alnum:]]", "_", reg)
  feat_name <- paste0("Reg_", clean_name, "1")
  
  if (feat_name %in% colnames(seu@meta.data)) {
    df <- data.frame(
      score = seu@meta.data[[feat_name]],
      Condition = seu$Condition,
      Cluster = seu$cl295v11SubShort
    ) %>% filter(!is.na(Condition))
    
    # Overall by condition
    p <- ggplot(df, aes(x = Condition, y = score, fill = Condition)) +
      geom_violin(alpha = 0.6, scale = "width") +
      geom_boxplot(width = 0.12, outlier.size = 0.2, alpha = 0.8) +
      scale_fill_manual(values = condition_colors) +
      labs(title = paste0(reg, " regulon"), y = "Module Score") +
      theme_cowplot() + NoLegend()
    
    violin_plots[[reg]] <- p
  }
}

if (length(violin_plots) > 0) {
  p_vln_all <- wrap_plots(violin_plots, ncol = 5)
  ggsave(file.path(result_dir, "regulon_violin_by_condition.pdf"), p_vln_all,
         width = 20, height = 4 * ceiling(length(violin_plots) / 5))
  ggsave(file.path(result_dir, "regulon_violin_by_condition.png"), p_vln_all,
         width = 20, height = 4 * ceiling(length(violin_plots) / 5), dpi = 200)
}

# --- F4. Condition-specific violin plots WITHIN macrophages only ---
cat("Generating macrophage-specific regulon violins by condition...\n")

seu_macro <- subset(seu, subset = cl295v11SubShort == "cM02" & !is.na(Condition))

macro_violin_plots <- list()
for (reg in key_regulons) {
  clean_name <- gsub("[^[:alnum:]]", "_", reg)
  feat_name <- paste0("Reg_", clean_name, "1")
  
  if (feat_name %in% colnames(seu_macro@meta.data)) {
    p <- VlnPlot(seu_macro, features = feat_name, group.by = "Condition",
                 pt.size = 0, cols = condition_colors) +
      labs(title = paste0(reg, " (Macrophages)")) +
      theme(plot.title = element_text(size = 10))
    macro_violin_plots[[reg]] <- p
  }
}

if (length(macro_violin_plots) > 0) {
  p_macro_vln <- wrap_plots(macro_violin_plots, ncol = 5)
  ggsave(file.path(result_dir, "regulon_violin_macrophage_by_condition.pdf"), p_macro_vln,
         width = 20, height = 4 * ceiling(length(macro_violin_plots) / 5))
  ggsave(file.path(result_dir, "regulon_violin_macrophage_by_condition.png"), p_macro_vln,
         width = 20, height = 4 * ceiling(length(macro_violin_plots) / 5), dpi = 200)
}

# --- F5. Statistical testing: regulon differences by condition ---
cat("\nStatistical testing: Regulon activity by condition in macrophages...\n")

stat_results <- data.frame()
for (reg in key_regulons) {
  clean_name <- gsub("[^[:alnum:]]", "_", reg)
  feat_name <- paste0("Reg_", clean_name, "1")
  
  if (feat_name %in% colnames(seu_macro@meta.data)) {
    scores <- seu_macro@meta.data[[feat_name]]
    conditions <- seu_macro$Condition
    
    # Kruskal-Wallis test
    kw <- kruskal.test(scores ~ conditions)
    
    # Pairwise Wilcoxon
    pw <- pairwise.wilcox.test(scores, conditions, p.adjust.method = "bonferroni")
    
    stat_results <- rbind(stat_results, data.frame(
      Regulon = reg,
      KW_pvalue = kw$p.value,
      N_vs_MMRp = pw$p.value["MMRp", "N"],
      N_vs_MMRd = pw$p.value["MMRd", "N"],
      MMRp_vs_MMRd = pw$p.value["MMRd", "MMRp"],
      Mean_N = mean(scores[conditions == "N"]),
      Mean_MMRp = mean(scores[conditions == "MMRp"]),
      Mean_MMRd = mean(scores[conditions == "MMRd"])
    ))
  }
}

stat_results <- stat_results %>% arrange(KW_pvalue)
write.csv(stat_results, file.path(result_dir, "regulon_stats_macrophage_by_condition.csv"), row.names = FALSE)

cat("\n--- Regulon Statistical Tests (Macrophages: N vs MMRp vs MMRd) ---\n")
print(stat_results)

# =========================================================================
# PART G: Regulon Specificity Score (RSS) + UMAP visualization
# =========================================================================
cat("\n=== Regulon Specificity & UMAP ===\n")

# RSS
compute_rss <- function(auc_mat, cluster_labels) {
  clusters <- unique(cluster_labels)
  rss_mat <- matrix(0, nrow = nrow(auc_mat), ncol = length(clusters),
                    dimnames = list(rownames(auc_mat), clusters))
  for (cl in clusters) {
    cl_indicator <- as.numeric(cluster_labels == cl) / sum(cluster_labels == cl)
    for (reg in rownames(auc_mat)) {
      reg_activity <- auc_mat[reg, ]
      reg_activity <- reg_activity / (sum(reg_activity) + 1e-10)
      m <- 0.5 * (reg_activity + cl_indicator)
      jsd <- 0.5 * sum(reg_activity * log2((reg_activity + 1e-10) / (m + 1e-10))) +
             0.5 * sum(cl_indicator * log2((cl_indicator + 1e-10) / (m + 1e-10)))
      rss_mat[reg, cl] <- 1 - sqrt(max(0, jsd))
    }
  }
  return(rss_mat)
}

cluster_labels <- meta_sub$cl295v11SubFull
rss <- compute_rss(regulonAUC, cluster_labels)
write.csv(rss, file.path(result_dir, "regulon_specificity_scores.csv"))

rss_long <- reshape2::melt(rss, varnames = c("Regulon", "Cluster"), value.name = "RSS")

p_rss <- ggplot(rss_long, aes(x = Cluster, y = Regulon, size = RSS, color = RSS)) +
  geom_point() +
  scale_color_viridis(option = "magma") +
  scale_size_continuous(range = c(0.5, 5)) +
  labs(title = "Regulon Specificity Score (RSS)") +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        axis.text.y = element_text(size = 8))

ggsave(file.path(result_dir, "regulon_specificity_dotplot.pdf"), p_rss, width = 10, height = 10)
ggsave(file.path(result_dir, "regulon_specificity_dotplot.png"), p_rss, width = 10, height = 10, dpi = 200)

# Regulon UMAP plots for key TFs
cat("Generating regulon UMAP plots...\n")
umap_reg_plots <- list()
for (reg in head(key_regulons, 6)) {
  clean_name <- gsub("[^[:alnum:]]", "_", reg)
  feat_name <- paste0("Reg_", clean_name, "1")
  if (feat_name %in% colnames(seu@meta.data)) {
    p <- FeaturePlot(seu, features = feat_name, pt.size = 0.2, raster = TRUE) +
      scale_color_viridis(option = "magma") +
      labs(title = paste0(reg, " regulon")) +
      theme_cowplot() + theme(plot.title = element_text(size = 10))
    umap_reg_plots[[reg]] <- p
  }
}

if (length(umap_reg_plots) > 0) {
  p_reg_umap <- wrap_plots(umap_reg_plots, ncol = 3)
  ggsave(file.path(result_dir, "regulon_umap_key_TFs.pdf"), p_reg_umap,
         width = 15, height = 5 * ceiling(length(umap_reg_plots) / 3))
  ggsave(file.path(result_dir, "regulon_umap_key_TFs.png"), p_reg_umap,
         width = 15, height = 5 * ceiling(length(umap_reg_plots) / 3), dpi = 200)
}

# --- Save ---
cat("\nSaving results...\n")
saveRDS(regulonAUC, file.path(result_dir, "regulon_AUC_matrix.rds"))
saveRDS(seu, file.path(result_dir, "seu_myeloid_regulons.rds"))

cat("\n Phase 8 complete!\n")
