# =============================================================================
# 05_differential_expression.R
# DE: Tumor vs Normal, MMRd vs MMRp in myeloid cells; macrophage subtype contrasts
# =============================================================================

library(Seurat)
library(ggplot2)
library(dplyr)
library(EnhancedVolcano)
library(ComplexHeatmap)
library(circlize)
library(patchwork)
library(cowplot)

set.seed(42)
cat("=== Phase 5: Differential Expression Analysis ===\n\n")

result_dir <- "results/05_de"
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

# --- 1. Load myeloid Seurat object ---
cat("Loading myeloid Seurat object...\n")
seu <- readRDS("results/04_myeloid/seu_myeloid.rds")

# Ensure RNA assay is default for DE
DefaultAssay(seu) <- "RNA"
seu <- NormalizeData(seu, verbose = FALSE)

# =========================================================================
# DE 1: Tumor vs Normal (all myeloid cells)
# =========================================================================
cat("\n--- DE: Tumor vs Normal (all myeloid) ---\n")

seu_tn <- subset(seu, subset = SPECIMEN_TYPE %in% c("T", "N"))
Idents(seu_tn) <- "SPECIMEN_TYPE"

de_tn <- FindMarkers(seu_tn, ident.1 = "T", ident.2 = "N",
                     test.use = "wilcox", logfc.threshold = 0.1,
                     min.pct = 0.1, verbose = TRUE)
de_tn$gene <- rownames(de_tn)

write.csv(de_tn, file.path(result_dir, "de_tumor_vs_normal_myeloid.csv"))
cat(sprintf("  Total DE genes: %d (up in Tumor: %d, up in Normal: %d)\n",
            nrow(de_tn),
            sum(de_tn$avg_log2FC > 0.25 & de_tn$p_val_adj < 0.05),
            sum(de_tn$avg_log2FC < -0.25 & de_tn$p_val_adj < 0.05)))

# Volcano plot
p_volcano_tn <- EnhancedVolcano(de_tn,
  lab = de_tn$gene, x = "avg_log2FC", y = "p_val_adj",
  title = "Myeloid: Tumor vs Normal",
  subtitle = "DE genes (Wilcoxon test)",
  pCutoff = 0.05, FCcutoff = 0.5,
  pointSize = 1.5, labSize = 3,
  col = c("grey30", "#2ECC71", "#3498DB", "#E74C3C"),
  drawConnectors = TRUE, widthConnectors = 0.3, maxoverlapsConnectors = 20)

ggsave(file.path(result_dir, "volcano_tumor_vs_normal.pdf"), p_volcano_tn, width = 10, height = 8)
ggsave(file.path(result_dir, "volcano_tumor_vs_normal.png"), p_volcano_tn, width = 10, height = 8, dpi = 200)

rm(seu_tn); gc()

# =========================================================================
# DE 2: MMRd vs MMRp (all myeloid cells)
# =========================================================================
cat("\n--- DE: MMRd vs MMRp (all myeloid) ---\n")

seu_mmr <- subset(seu, subset = MMRStatus %in% c("MMRd", "MMRp"))
Idents(seu_mmr) <- "MMRStatus"

de_mmr <- FindMarkers(seu_mmr, ident.1 = "MMRd", ident.2 = "MMRp",
                      test.use = "wilcox", logfc.threshold = 0.1,
                      min.pct = 0.1, verbose = TRUE)
de_mmr$gene <- rownames(de_mmr)

write.csv(de_mmr, file.path(result_dir, "de_mmrd_vs_mmrp_myeloid.csv"))
cat(sprintf("  Total DE genes: %d (up in MMRd: %d, up in MMRp: %d)\n",
            nrow(de_mmr),
            sum(de_mmr$avg_log2FC > 0.25 & de_mmr$p_val_adj < 0.05),
            sum(de_mmr$avg_log2FC < -0.25 & de_mmr$p_val_adj < 0.05)))

p_volcano_mmr <- EnhancedVolcano(de_mmr,
  lab = de_mmr$gene, x = "avg_log2FC", y = "p_val_adj",
  title = "Myeloid: MMRd vs MMRp",
  subtitle = "DE genes (Wilcoxon test)",
  pCutoff = 0.05, FCcutoff = 0.5,
  pointSize = 1.5, labSize = 3,
  col = c("grey30", "#2ECC71", "#3498DB", "#E74C3C"),
  drawConnectors = TRUE, widthConnectors = 0.3, maxoverlapsConnectors = 20)

ggsave(file.path(result_dir, "volcano_mmrd_vs_mmrp.pdf"), p_volcano_mmr, width = 10, height = 8)
ggsave(file.path(result_dir, "volcano_mmrd_vs_mmrp.png"), p_volcano_mmr, width = 10, height = 8, dpi = 200)

rm(seu_mmr); gc()

# =========================================================================
# DE 3: Per-subcluster Tumor vs Normal
# =========================================================================
cat("\n--- DE: Tumor vs Normal per myeloid subcluster ---\n")

subclusters <- unique(seu$cl295v11SubShort)
de_per_cluster <- list()

for (sc in subclusters) {
  cat(sprintf("  Processing %s...\n", sc))
  seu_sub <- subset(seu, subset = cl295v11SubShort == sc & SPECIMEN_TYPE %in% c("T", "N"))
  
  n_t <- sum(seu_sub$SPECIMEN_TYPE == "T")
  n_n <- sum(seu_sub$SPECIMEN_TYPE == "N")
  
  if (n_t >= 20 & n_n >= 20) {
    Idents(seu_sub) <- "SPECIMEN_TYPE"
    de <- tryCatch({
      FindMarkers(seu_sub, ident.1 = "T", ident.2 = "N",
                  test.use = "wilcox", logfc.threshold = 0.1,
                  min.pct = 0.1, verbose = FALSE)
    }, error = function(e) NULL)
    
    if (!is.null(de)) {
      de$gene <- rownames(de)
      de$cluster <- sc
      de_per_cluster[[sc]] <- de
      write.csv(de, file.path(result_dir, paste0("de_TvsN_", sc, ".csv")))
      cat(sprintf("    %d DE genes found\n", nrow(de)))
    }
  } else {
    cat(sprintf("    Skipped (T=%d, N=%d cells)\n", n_t, n_n))
  }
}

# Combine all per-cluster DE results
if (length(de_per_cluster) > 0) {
  de_all <- bind_rows(de_per_cluster)
  write.csv(de_all, file.path(result_dir, "de_TvsN_all_subclusters.csv"), row.names = FALSE)
}

# =========================================================================
# DE 4: Macrophage subtype contrasts
# =========================================================================
cat("\n--- DE: Macrophage subtype contrasts ---\n")

# cM02 (Macrophage-like) vs cM01 (Monocyte)
cat("  cM02 (Macrophage) vs cM01 (Monocyte)...\n")
seu_mac_mono <- subset(seu, subset = cl295v11SubShort %in% c("cM02", "cM01"))
Idents(seu_mac_mono) <- "cl295v11SubShort"

de_mac_mono <- FindMarkers(seu_mac_mono, ident.1 = "cM02", ident.2 = "cM01",
                           test.use = "wilcox", logfc.threshold = 0.1,
                           min.pct = 0.1, verbose = FALSE)
de_mac_mono$gene <- rownames(de_mac_mono)
write.csv(de_mac_mono, file.path(result_dir, "de_macrophage_vs_monocyte.csv"))

p_volcano_mac <- EnhancedVolcano(de_mac_mono,
  lab = de_mac_mono$gene, x = "avg_log2FC", y = "p_val_adj",
  title = "Macrophage (cM02) vs Monocyte (cM01)",
  pCutoff = 0.05, FCcutoff = 0.5,
  pointSize = 1.5, labSize = 3,
  drawConnectors = TRUE, widthConnectors = 0.3)

ggsave(file.path(result_dir, "volcano_macrophage_vs_monocyte.pdf"), p_volcano_mac, width = 10, height = 8)
ggsave(file.path(result_dir, "volcano_macrophage_vs_monocyte.png"), p_volcano_mac, width = 10, height = 8, dpi = 200)

# =========================================================================
# DE 5: Top DE genes heatmap
# =========================================================================
cat("\nGenerating DE heatmap...\n")

# Top DE genes from Tumor vs Normal
top_up   <- de_tn %>% filter(p_val_adj < 0.05) %>% slice_max(n = 30, order_by = avg_log2FC)
top_down <- de_tn %>% filter(p_val_adj < 0.05) %>% slice_min(n = 30, order_by = avg_log2FC)
top_genes <- c(top_up$gene, top_down$gene)

avg_expr <- AverageExpression(seu, features = top_genes,
                               group.by = "cl295v11SubFull")
avg_mat <- as.matrix(avg_expr$RNA)
avg_mat_scaled <- t(scale(t(avg_mat)))

pdf(file.path(result_dir, "de_heatmap_TvsN_top_genes.pdf"), width = 10, height = 12)
ht <- Heatmap(avg_mat_scaled, name = "Z-score",
              col = colorRamp2(c(-2, 0, 2), c("#3498DB", "white", "#E74C3C")),
              row_names_gp = gpar(fontsize = 7),
              column_names_gp = gpar(fontsize = 8),
              column_title = "Top DE Genes: Tumor vs Normal (Myeloid)")
draw(ht)
dev.off()

cat("\n Phase 5 complete!\n")
