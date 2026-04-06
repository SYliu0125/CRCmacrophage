# =============================================================================
# 09_comparison_published.R
# Compare reanalysis results with published findings from Pelka et al. 2021
# =============================================================================

library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(cowplot)
library(ComplexHeatmap)
library(circlize)
library(viridis)

set.seed(42)
cat("=== Phase 9: Comparison to Published Findings ===\n\n")

result_dir <- "results/09_comparison"
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

# --- 1. Load objects ---
cat("Loading Seurat objects...\n")
seu_full <- readRDS("results/02_processing/seu_full_processed.rds")
seu_myeloid <- readRDS("results/04_myeloid/seu_myeloid.rds")

# =========================================================================
# 1. Reproduce key published cell type proportions
# =========================================================================
cat("\n--- Comparison 1: Cell Type Proportions ---\n")

# Published approximate proportions from the paper (Extended Data Fig 1):
published_proportions <- data.frame(
  CellType = c("Epi", "TNKILC", "Myeloid", "Plasma", "B", "Strom", "Mast"),
  Published_pct = c(45.5, 20.8, 11.4, 10.2, 6.9, 4.2, 1.0)
)

# Our proportions
our_counts <- table(seu_full$clTopLevel)
our_proportions <- data.frame(
  CellType = names(our_counts),
  Our_pct = as.numeric(100 * our_counts / sum(our_counts))
)

comparison <- merge(published_proportions, our_proportions, by = "CellType")
comparison$Difference <- comparison$Our_pct - comparison$Published_pct

cat("Cell type proportion comparison:\n")
print(comparison)
write.csv(comparison, file.path(result_dir, "celltype_proportion_comparison.csv"), row.names = FALSE)

# Plot comparison
comp_long <- comparison %>%
  select(CellType, Published_pct, Our_pct) %>%
  pivot_longer(cols = c(Published_pct, Our_pct),
               names_to = "Source", values_to = "Percentage") %>%
  mutate(Source = ifelse(Source == "Published_pct", "Published", "Our Analysis"))

p1 <- ggplot(comp_long, aes(x = CellType, y = Percentage, fill = Source)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("Published" = "#3498DB", "Our Analysis" = "#E74C3C")) +
  labs(title = "Cell Type Proportions: Published vs Our Analysis",
       x = "Cell Type", y = "Percentage (%)") +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave(file.path(result_dir, "proportion_comparison.pdf"), p1, width = 8, height = 5)
ggsave(file.path(result_dir, "proportion_comparison.png"), p1, width = 8, height = 5, dpi = 200)

# =========================================================================
# 2. Validate myeloid subtype markers (Published Table)
# =========================================================================
cat("\n--- Comparison 2: Myeloid Subtype Marker Validation ---\n")

# Published key markers per myeloid subcluster from Pelka et al. 2021
published_markers <- list(
  "cM01 (Monocyte)" = c("S100A8", "S100A9", "VCAN", "FCN1", "CD14", "LYZ"),
  "cM02 (Macrophage-like)" = c("C1QA", "C1QB", "C1QC", "APOE", "CD68", "CD163", "MRC1"),
  "cM03 (DC1)" = c("CLEC9A", "XCR1", "CADM1", "BATF3", "IDO1"),
  "cM04 (DC2)" = c("CD1C", "FCER1A", "CLEC10A", "CD1E"),
  "cM05 (DC2 C1Q+)" = c("C1QA", "C1QB", "CD1C", "FCER1A"),
  "cM06 (DC IL22RA2)" = c("LAMP3", "CCR7", "FSCN1", "IL22RA2"),
  "cM07 (pDC)" = c("LILRA4", "IRF7", "TCF4", "ITM2C", "PLD4"),
  "cM09 (mregDC)" = c("CCR7", "LAMP3", "CCL19", "FSCN1"),
  "cM10 (Granulocyte)" = c("CSF3R", "CXCR2", "S100A8", "S100A9")
)

# Check expression of published markers in our data
DefaultAssay(seu_myeloid) <- "RNA"
Idents(seu_myeloid) <- "cl295v11SubFull"

# Dot plot of published markers
all_published_markers <- unique(unlist(published_markers))
all_published_markers <- all_published_markers[all_published_markers %in% rownames(seu_myeloid)]

p2 <- DotPlot(seu_myeloid, features = all_published_markers, group.by = "cl295v11SubFull") +
  coord_flip() +
  scale_color_viridis(option = "magma") +
  labs(title = "Published Myeloid Markers in Our Data") +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        axis.text.y = element_text(size = 7))

ggsave(file.path(result_dir, "published_markers_validation.pdf"), p2, width = 10, height = 12)
ggsave(file.path(result_dir, "published_markers_validation.png"), p2, width = 10, height = 12, dpi = 200)

# --- Quantitative marker validation ---
cat("\nQuantitative marker validation:\n")

# For each cluster, check if published markers are significantly enriched
marker_validation <- data.frame()

for (cluster in names(published_markers)) {
  markers <- published_markers[[cluster]]
  markers <- markers[markers %in% rownames(seu_myeloid)]
  
  if (length(markers) > 0) {
    # Get average expression in target cluster vs others
    avg_expr <- AverageExpression(seu_myeloid, features = markers, group.by = "cl295v11SubFull")
    avg_mat <- as.matrix(avg_expr$RNA)
    
    if (cluster %in% colnames(avg_mat)) {
      target_expr <- mean(avg_mat[, cluster])
      other_expr  <- mean(avg_mat[, colnames(avg_mat) != cluster])
      enrichment  <- target_expr / (other_expr + 0.01)
      
      marker_validation <- rbind(marker_validation, data.frame(
        Cluster = cluster,
        N_markers = length(markers),
        Avg_expr_target = round(target_expr, 3),
        Avg_expr_other = round(other_expr, 3),
        Enrichment_fold = round(enrichment, 2),
        Validated = ifelse(enrichment > 1.5, "YES", "PARTIAL")
      ))
      
      cat(sprintf("  %-30s markers=%d, enrichment=%.1fx -> %s\n",
                  cluster, length(markers), enrichment,
                  ifelse(enrichment > 1.5, "VALIDATED", "PARTIAL")))
    }
  }
}

write.csv(marker_validation, file.path(result_dir, "marker_validation_results.csv"), row.names = FALSE)

# =========================================================================
# 3. Compare our UMAP with published cluster positions
# =========================================================================
cat("\n--- Comparison 3: UMAP Cluster Separation ---\n")

# Silhouette-like metric: measure how well clusters separate
cat("Assessing cluster separation...\n")

umap_coords <- Embeddings(seu_myeloid, "umap")
clusters <- seu_myeloid$cl295v11SubShort

# Compute cluster centroids and inter-cluster distances
centroids <- aggregate(umap_coords, by = list(Cluster = clusters), FUN = mean)
rownames(centroids) <- centroids$Cluster
centroids$Cluster <- NULL

centroid_dist <- as.matrix(dist(centroids))
cat("  Centroid distance matrix:\n")
print(round(centroid_dist, 2))

write.csv(centroid_dist, file.path(result_dir, "cluster_centroid_distances.csv"))

# =========================================================================
# 4. Compare DE gene overlap with published
# =========================================================================
cat("\n--- Comparison 4: DE Gene Concordance ---\n")

# Load our DE results
if (file.exists("results/05_de/de_tumor_vs_normal_myeloid.csv")) {
  de_tn <- read.csv("results/05_de/de_tumor_vs_normal_myeloid.csv", row.names = 1)
  
  # Key published upregulated genes in tumor myeloid cells
  published_tumor_up <- c("SPP1", "VEGFA", "MARCO", "TREM2", "GPNMB",
                           "FN1", "MMP9", "IL1B", "CXCL8", "CCL2",
                           "HIF1A", "LDHA", "SLC2A1")
  
  published_tumor_up_in_data <- published_tumor_up[published_tumor_up %in% de_tn$gene]
  
  our_sig_up <- de_tn %>% 
    filter(p_val_adj < 0.05, avg_log2FC > 0.25) %>%
    arrange(desc(avg_log2FC))
  
  overlap <- intersect(published_tumor_up_in_data, our_sig_up$gene)
  
  cat(sprintf("  Published tumor-upregulated genes: %d\n", length(published_tumor_up)))
  cat(sprintf("  In our data: %d\n", length(published_tumor_up_in_data)))
  cat(sprintf("  Also significant in our DE: %d (%.0f%%)\n",
              length(overlap), 100 * length(overlap) / length(published_tumor_up_in_data)))
  cat(sprintf("  Concordant genes: %s\n", paste(overlap, collapse = ", ")))
  
  # Check direction concordance
  concordance_df <- data.frame(
    Gene = published_tumor_up_in_data,
    Published_direction = "UP in Tumor"
  )
  
  for (gene in published_tumor_up_in_data) {
    if (gene %in% rownames(de_tn)) {
      concordance_df$Our_log2FC[concordance_df$Gene == gene] <- de_tn[gene, "avg_log2FC"]
      concordance_df$Our_padj[concordance_df$Gene == gene] <- de_tn[gene, "p_val_adj"]
      concordance_df$Our_direction[concordance_df$Gene == gene] <- 
        ifelse(de_tn[gene, "avg_log2FC"] > 0, "UP", "DOWN")
      concordance_df$Concordant[concordance_df$Gene == gene] <- 
        de_tn[gene, "avg_log2FC"] > 0
    }
  }
  
  write.csv(concordance_df, file.path(result_dir, "de_concordance_with_published.csv"), row.names = FALSE)
  cat("\nDE concordance table:\n")
  print(concordance_df)
}

# =========================================================================
# 5. Summary report
# =========================================================================
cat("\n\n========================================\n")
cat("     COMPARISON SUMMARY\n")
cat("========================================\n\n")

cat("1. CELL TYPE PROPORTIONS:\n")
cat(sprintf("   Mean absolute difference: %.1f%%\n", mean(abs(comparison$Difference))))
cat(sprintf("   Max difference: %.1f%% (%s)\n",
            max(abs(comparison$Difference)), 
            comparison$CellType[which.max(abs(comparison$Difference))]))

cat("\n2. MYELOID MARKER VALIDATION:\n")
cat(sprintf("   Clusters validated: %d / %d\n",
            sum(marker_validation$Validated == "YES"),
            nrow(marker_validation)))

cat("\n3. DE GENE CONCORDANCE:\n")
if (exists("concordance_df")) {
  cat(sprintf("   Direction agreement: %d / %d (%.0f%%)\n",
              sum(concordance_df$Concordant, na.rm = TRUE),
              nrow(concordance_df),
              100 * sum(concordance_df$Concordant, na.rm = TRUE) / nrow(concordance_df)))
}

cat("\n Phase 9 complete!\n")
