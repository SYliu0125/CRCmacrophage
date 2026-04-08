# =============================================================================
# 05_differential_expression.R
# DE: Tumor vs Normal, MMRd vs MMRp in myeloid cells; macrophage subtype contrasts
#
# METHOD NOTE:
#   Condition comparisons (T vs N, MMRd vs MMRp) use PSEUDOBULK DESeq2:
#     counts are aggregated per patient (PID) × condition, then DESeq2 is run
#     on those ~30-90 patient-level samples. This avoids pseudoreplication
#     (cells from the same patient are not independent observations).
#   Cluster marker contrasts (cM02 vs cM01) retain Wilcoxon via FindMarkers
#     because these compare cell types within the same sample, not conditions
#     across patients.
# =============================================================================

library(Seurat)
library(DESeq2)
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
DefaultAssay(seu) <- "RNA"
seu <- NormalizeData(seu, verbose = FALSE)

# =========================================================================
# PSEUDOBULK HELPER FUNCTIONS
# =========================================================================

# Aggregate raw counts per PID x condition into a pseudobulk matrix.
# Returns list(counts = matrix, coldata = data.frame)
make_pseudobulk <- function(seu_sub, condition_col, min_cells = 10) {
  counts <- GetAssayData(seu_sub, slot = "counts", assay = "RNA")
  meta   <- seu_sub@meta.data[, c("PID", condition_col), drop = FALSE]
  meta$pb_sample <- paste(meta$PID, meta[[condition_col]], sep = "__")

  groups   <- unique(meta$pb_sample)
  pb_list  <- lapply(groups, function(g) {
    cells <- rownames(meta)[meta$pb_sample == g]
    if (length(cells) < min_cells) return(NULL)
    Matrix::rowSums(counts[, cells, drop = FALSE])
  })
  names(pb_list) <- groups
  pb_list <- pb_list[!sapply(pb_list, is.null)]

  pb_mat  <- do.call(cbind, pb_list)

  coldata <- data.frame(
    sample    = names(pb_list),
    PID       = sub("__.*$",  "", names(pb_list)),
    condition = sub("^.*__",  "", names(pb_list)),
    row.names = names(pb_list),
    stringsAsFactors = FALSE
  )
  colnames(coldata)[colnames(coldata) == "condition"] <- condition_col

  list(counts = as.matrix(pb_mat), coldata = coldata)
}

# Run DESeq2 on a pseudobulk object.
# level1 = numerator (e.g. "T"), level2 = denominator (e.g. "N")
run_deseq2 <- function(pb, condition_col, level1, level2) {
  coldata <- pb$coldata
  coldata[[condition_col]] <- factor(coldata[[condition_col]],
                                     levels = c(level2, level1))

  # Keep only samples belonging to the two levels being compared
  keep_samples <- coldata[[condition_col]] %in% c(level1, level2)
  coldata  <- coldata[keep_samples, , drop = FALSE]
  pb_mat   <- pb$counts[, rownames(coldata), drop = FALSE]

  # Require at least 3 samples per group
  n1 <- sum(coldata[[condition_col]] == level1)
  n2 <- sum(coldata[[condition_col]] == level2)
  if (n1 < 3 | n2 < 3) {
    warning(sprintf("Too few pseudobulk samples: %s n=%d, %s n=%d. Skipping.",
                    level1, n1, level2, n2))
    return(NULL)
  }

  # Filter lowly expressed genes (expressed in >= 2 samples)
  keep_genes <- rowSums(pb_mat >= 1) >= 2
  pb_mat <- pb_mat[keep_genes, ]

  design_formula <- as.formula(paste0("~", condition_col))
  dds <- DESeqDataSetFromMatrix(countData = pb_mat,
                                colData   = coldata,
                                design    = design_formula)
  dds <- DESeq(dds, quiet = TRUE)
  res <- results(dds, contrast = c(condition_col, level1, level2))

  df <- as.data.frame(res)
  df$gene <- rownames(df)
  df
}

# =========================================================================
# DE 1: Tumor vs Normal (all myeloid cells) — PSEUDOBULK
# =========================================================================
cat("\n--- DE: Tumor vs Normal (all myeloid) — Pseudobulk DESeq2 ---\n")

seu_tn  <- subset(seu, subset = SPECIMEN_TYPE %in% c("T", "N"))
pb_tn   <- make_pseudobulk(seu_tn, condition_col = "SPECIMEN_TYPE")

cat(sprintf("  Pseudobulk samples: %d total (T=%d, N=%d)\n",
            nrow(pb_tn$coldata),
            sum(pb_tn$coldata$SPECIMEN_TYPE == "T"),
            sum(pb_tn$coldata$SPECIMEN_TYPE == "N")))

de_tn <- run_deseq2(pb_tn, condition_col = "SPECIMEN_TYPE",
                    level1 = "T", level2 = "N")

write.csv(de_tn, file.path(result_dir, "de_tumor_vs_normal_myeloid.csv"), row.names = FALSE)
cat(sprintf("  DE genes (|log2FC| > 0.5, padj < 0.05): up in Tumor=%d, up in Normal=%d\n",
            sum(de_tn$log2FoldChange >  0.5 & de_tn$padj < 0.05, na.rm = TRUE),
            sum(de_tn$log2FoldChange < -0.5 & de_tn$padj < 0.05, na.rm = TRUE)))

p_volcano_tn <- EnhancedVolcano(de_tn,
  lab = de_tn$gene, x = "log2FoldChange", y = "padj",
  title    = "Myeloid: Tumor vs Normal",
  subtitle = "Pseudobulk DESeq2",
  pCutoff  = 0.05, FCcutoff = 0.5,
  pointSize = 1.5, labSize = 3,
  col = c("grey30", "#2ECC71", "#3498DB", "#E74C3C"),
  drawConnectors = TRUE, widthConnectors = 0.3, maxoverlapsConnectors = 20)

ggsave(file.path(result_dir, "volcano_tumor_vs_normal.pdf"), p_volcano_tn, width = 10, height = 8)
ggsave(file.path(result_dir, "volcano_tumor_vs_normal.png"), p_volcano_tn, width = 10, height = 8, dpi = 200)

rm(seu_tn); gc()

# =========================================================================
# DE 2: MMRd vs MMRp (all myeloid cells) — PSEUDOBULK
# =========================================================================
cat("\n--- DE: MMRd vs MMRp (all myeloid) — Pseudobulk DESeq2 ---\n")

seu_mmr <- subset(seu, subset = MMRStatus %in% c("MMRd", "MMRp"))
pb_mmr  <- make_pseudobulk(seu_mmr, condition_col = "MMRStatus")

cat(sprintf("  Pseudobulk samples: %d total (MMRd=%d, MMRp=%d)\n",
            nrow(pb_mmr$coldata),
            sum(pb_mmr$coldata$MMRStatus == "MMRd"),
            sum(pb_mmr$coldata$MMRStatus == "MMRp")))

de_mmr <- run_deseq2(pb_mmr, condition_col = "MMRStatus",
                     level1 = "MMRd", level2 = "MMRp")

write.csv(de_mmr, file.path(result_dir, "de_mmrd_vs_mmrp_myeloid.csv"), row.names = FALSE)
cat(sprintf("  DE genes (|log2FC| > 0.5, padj < 0.05): up in MMRd=%d, up in MMRp=%d\n",
            sum(de_mmr$log2FoldChange >  0.5 & de_mmr$padj < 0.05, na.rm = TRUE),
            sum(de_mmr$log2FoldChange < -0.5 & de_mmr$padj < 0.05, na.rm = TRUE)))

p_volcano_mmr <- EnhancedVolcano(de_mmr,
  lab = de_mmr$gene, x = "log2FoldChange", y = "padj",
  title    = "Myeloid: MMRd vs MMRp",
  subtitle = "Pseudobulk DESeq2",
  pCutoff  = 0.05, FCcutoff = 0.5,
  pointSize = 1.5, labSize = 3,
  col = c("grey30", "#2ECC71", "#3498DB", "#E74C3C"),
  drawConnectors = TRUE, widthConnectors = 0.3, maxoverlapsConnectors = 20)

ggsave(file.path(result_dir, "volcano_mmrd_vs_mmrp.pdf"), p_volcano_mmr, width = 10, height = 8)
ggsave(file.path(result_dir, "volcano_mmrd_vs_mmrp.png"), p_volcano_mmr, width = 10, height = 8, dpi = 200)

rm(seu_mmr); gc()

# =========================================================================
# DE 3: Per-subcluster Tumor vs Normal — PSEUDOBULK
# =========================================================================
cat("\n--- DE: Tumor vs Normal per myeloid subcluster — Pseudobulk DESeq2 ---\n")

subclusters      <- unique(seu$cl295v11SubShort)
de_per_cluster   <- list()

for (sc in subclusters) {
  cat(sprintf("  Processing %s...\n", sc))
  seu_sub <- subset(seu, subset = cl295v11SubShort == sc &
                                  SPECIMEN_TYPE %in% c("T", "N"))

  pb_sub <- make_pseudobulk(seu_sub, condition_col = "SPECIMEN_TYPE")

  n_t <- sum(pb_sub$coldata$SPECIMEN_TYPE == "T")
  n_n <- sum(pb_sub$coldata$SPECIMEN_TYPE == "N")

  if (n_t >= 3 & n_n >= 3) {
    de <- tryCatch(
      run_deseq2(pb_sub, condition_col = "SPECIMEN_TYPE",
                 level1 = "T", level2 = "N"),
      error = function(e) { message("    Error: ", e$message); NULL }
    )
    if (!is.null(de)) {
      de$cluster <- sc
      de_per_cluster[[sc]] <- de
      write.csv(de, file.path(result_dir, paste0("de_TvsN_", sc, ".csv")),
                row.names = FALSE)
      cat(sprintf("    %d DE genes found (padj < 0.05)\n",
                  sum(de$padj < 0.05, na.rm = TRUE)))
    }
  } else {
    cat(sprintf("    Skipped (pseudobulk T=%d, N=%d patients)\n", n_t, n_n))
  }
}

if (length(de_per_cluster) > 0) {
  de_all <- bind_rows(de_per_cluster)
  write.csv(de_all, file.path(result_dir, "de_TvsN_all_subclusters.csv"),
            row.names = FALSE)
}

# =========================================================================
# DE 4: Macrophage subtype contrasts — Wilcoxon (cell-type comparison, not
#        condition comparison, so pseudobulk is not required here)
# =========================================================================
cat("\n--- DE: Macrophage subtype contrasts (Wilcoxon) ---\n")

seu_mac_mono <- subset(seu, subset = cl295v11SubShort %in% c("cM02", "cM01"))
Idents(seu_mac_mono) <- "cl295v11SubShort"

de_mac_mono <- FindMarkers(seu_mac_mono, ident.1 = "cM02", ident.2 = "cM01",
                           test.use = "wilcox", logfc.threshold = 0.1,
                           min.pct = 0.1, verbose = FALSE)
de_mac_mono$gene <- rownames(de_mac_mono)
write.csv(de_mac_mono, file.path(result_dir, "de_macrophage_vs_monocyte.csv"))

p_volcano_mac <- EnhancedVolcano(de_mac_mono,
  lab = de_mac_mono$gene, x = "avg_log2FC", y = "p_val_adj",
  title    = "Macrophage (cM02) vs Monocyte (cM01)",
  subtitle = "Wilcoxon (cell-type marker comparison)",
  pCutoff  = 0.05, FCcutoff = 0.5,
  pointSize = 1.5, labSize = 3,
  drawConnectors = TRUE, widthConnectors = 0.3)

ggsave(file.path(result_dir, "volcano_macrophage_vs_monocyte.pdf"), p_volcano_mac, width = 10, height = 8)
ggsave(file.path(result_dir, "volcano_macrophage_vs_monocyte.png"), p_volcano_mac, width = 10, height = 8, dpi = 200)

# =========================================================================
# DE 5: Top DE genes heatmap (Tumor vs Normal)
# =========================================================================
cat("\nGenerating DE heatmap...\n")

top_up   <- de_tn %>% filter(padj < 0.05) %>% slice_max(n = 30, order_by = log2FoldChange)
top_down <- de_tn %>% filter(padj < 0.05) %>% slice_min(n = 30, order_by = log2FoldChange)
top_genes <- c(top_up$gene, top_down$gene)
top_genes <- top_genes[top_genes %in% rownames(seu)]

avg_expr       <- AverageExpression(seu, features = top_genes,
                                    group.by = "cl295v11SubFull")
avg_mat        <- as.matrix(avg_expr$RNA)
avg_mat_scaled <- t(scale(t(avg_mat)))

pdf(file.path(result_dir, "de_heatmap_TvsN_top_genes.pdf"), width = 10, height = 12)
ht <- Heatmap(avg_mat_scaled, name = "Z-score",
              col = colorRamp2(c(-2, 0, 2), c("#3498DB", "white", "#E74C3C")),
              row_names_gp = gpar(fontsize = 7),
              column_names_gp = gpar(fontsize = 8),
              column_title = "Top DE Genes: Tumor vs Normal (Pseudobulk DESeq2)")
draw(ht)
dev.off()

cat("\n Phase 5 complete!\n")
