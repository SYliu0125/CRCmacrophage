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
#
# ANALYSIS STEPS:
#   Step 1 — Tumor vs Normal, all myeloid             → pseudobulk DESeq2
#   Step 2 — MMRd vs MMRp, all myeloid                → pseudobulk DESeq2
#   Step 3 — Tumor vs Normal, per myeloid subcluster  → pseudobulk DESeq2
#   Step 4 — Macrophage (cM02) vs Monocyte (cM01)     → Wilcoxon (FindMarkers)
#   Step 5 — Top DE genes heatmap                     → uses Step 1 output
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

# =============================================================================
# LOAD DATA
# Input:  results/04_myeloid/seu_myeloid.rds
# Output: seu  — myeloid Seurat object with RNA assay normalized
# =============================================================================
cat("Loading myeloid Seurat object...\n")
seu <- readRDS("results/04_myeloid/seu_myeloid.rds")
DefaultAssay(seu) <- "RNA"
seu <- NormalizeData(seu, verbose = FALSE)


# =============================================================================
# STEP 1: Tumor vs Normal — all myeloid cells (pseudobulk DESeq2)
# Input:  seu,  SPECIMEN_TYPE column  ("T" = tumor, "N" = normal)
# Output: de_tumor_vs_normal_myeloid.csv
#         volcano_tumor_vs_normal.pdf / .png
#         de_tn  — kept in memory, used again in Step 5
# =============================================================================
cat("\n--- STEP 1: Tumor vs Normal (all myeloid) ---\n")

# 1a. Keep only tumor (T) and normal (N) cells
seu_tn <- subset(seu, subset = SPECIMEN_TYPE %in% c("T", "N"))

# 1b. Extract the raw (integer) count matrix: genes × cells
counts_tn <- GetAssayData(seu_tn, slot = "counts", assay = "RNA")

# 1c. Create a pseudobulk group label "PID__condition" for each cell.
#     Cells from the same patient × same condition will be summed together.
meta_tn <- seu_tn@meta.data[, c("PID", "SPECIMEN_TYPE")]
meta_tn$pb_group <- paste(meta_tn$PID, meta_tn$SPECIMEN_TYPE, sep = "__")

# 1d. For each unique PID×condition group, sum raw counts across its cells.
#     Groups with fewer than 10 cells are dropped (too sparse to be reliable).
groups_tn <- unique(meta_tn$pb_group)
pb_list_tn <- lapply(groups_tn, function(g) {
  cells_in_group <- rownames(meta_tn)[meta_tn$pb_group == g]
  if (length(cells_in_group) < 10) return(NULL)
  Matrix::rowSums(counts_tn[, cells_in_group, drop = FALSE])
})
names(pb_list_tn) <- groups_tn
pb_list_tn <- pb_list_tn[!sapply(pb_list_tn, is.null)]  # remove dropped groups

# 1e. Assemble into a gene × pseudobulk-sample count matrix
pb_counts_tn <- do.call(cbind, pb_list_tn)

# 1f. Build a sample metadata table (one row per pseudobulk sample)
#     DESeq2 needs this to know which condition each column belongs to.
pb_meta_tn <- data.frame(
  sample         = names(pb_list_tn),
  PID            = sub("__.*$", "", names(pb_list_tn)),   # extract patient ID
  SPECIMEN_TYPE  = sub("^.*__", "", names(pb_list_tn)),   # extract T or N
  row.names      = names(pb_list_tn),
  stringsAsFactors = FALSE
)
pb_meta_tn$SPECIMEN_TYPE <- factor(pb_meta_tn$SPECIMEN_TYPE, levels = c("N", "T"))
# N is the reference level (denominator); T is the numerator

cat(sprintf("  Pseudobulk samples: %d total (T=%d, N=%d)\n",
            nrow(pb_meta_tn),
            sum(pb_meta_tn$SPECIMEN_TYPE == "T"),
            sum(pb_meta_tn$SPECIMEN_TYPE == "N")))

# 1g. Require at least 3 pseudobulk samples per condition before running DESeq2
stopifnot(sum(pb_meta_tn$SPECIMEN_TYPE == "T") >= 3,
          sum(pb_meta_tn$SPECIMEN_TYPE == "N") >= 3)

# 1h. Filter out genes expressed in fewer than 2 pseudobulk samples (too sparse)
keep_genes_tn <- rowSums(pb_counts_tn >= 1) >= 2
pb_counts_tn  <- pb_counts_tn[keep_genes_tn, ]

# 1i. Create DESeq2 dataset object
dds_tn <- DESeqDataSetFromMatrix(
  countData = as.matrix(pb_counts_tn),
  colData   = pb_meta_tn,
  design    = ~ SPECIMEN_TYPE    # model: expression ~ condition
)

# 1j. Run DESeq2: estimate size factors, dispersion, and fit the model
dds_tn <- DESeq(dds_tn, quiet = TRUE)

# 1k. Extract results: log2FC and adjusted p-value for T vs N
de_tn       <- as.data.frame(results(dds_tn, contrast = c("SPECIMEN_TYPE", "T", "N")))
de_tn$gene  <- rownames(de_tn)

write.csv(de_tn, file.path(result_dir, "de_tumor_vs_normal_myeloid.csv"), row.names = FALSE)
cat(sprintf("  DE genes (|log2FC| > 0.5, padj < 0.05): up in Tumor=%d, up in Normal=%d\n",
            sum(de_tn$log2FoldChange >  0.5 & de_tn$padj < 0.05, na.rm = TRUE),
            sum(de_tn$log2FoldChange < -0.5 & de_tn$padj < 0.05, na.rm = TRUE)))

# 1l. Volcano plot
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

rm(seu_tn, counts_tn, meta_tn, pb_list_tn, pb_counts_tn, pb_meta_tn, dds_tn); gc()


# =============================================================================
# STEP 2: MMRd vs MMRp — all myeloid cells (pseudobulk DESeq2)
# Input:  seu,  MMRStatus column  ("MMRd" = mismatch repair deficient,
#                                  "MMRp" = mismatch repair proficient)
# Output: de_mmrd_vs_mmrp_myeloid.csv
#         volcano_mmrd_vs_mmrp.pdf / .png
# =============================================================================
cat("\n--- STEP 2: MMRd vs MMRp (all myeloid) ---\n")

# 2a. Keep only cells with known MMR status (excludes normal tissue)
seu_mmr <- subset(seu, subset = MMRStatus %in% c("MMRd", "MMRp"))

# 2b. Extract raw count matrix
counts_mmr <- GetAssayData(seu_mmr, slot = "counts", assay = "RNA")

# 2c. Pseudobulk group label: PID × MMRStatus
meta_mmr <- seu_mmr@meta.data[, c("PID", "MMRStatus")]
meta_mmr$pb_group <- paste(meta_mmr$PID, meta_mmr$MMRStatus, sep = "__")

# 2d. Sum counts per PID × MMRStatus group (drop groups with < 10 cells)
groups_mmr <- unique(meta_mmr$pb_group)
pb_list_mmr <- lapply(groups_mmr, function(g) {
  cells_in_group <- rownames(meta_mmr)[meta_mmr$pb_group == g]
  if (length(cells_in_group) < 10) return(NULL)
  Matrix::rowSums(counts_mmr[, cells_in_group, drop = FALSE])
})
names(pb_list_mmr) <- groups_mmr
pb_list_mmr <- pb_list_mmr[!sapply(pb_list_mmr, is.null)]

# 2e. Assemble pseudobulk count matrix
pb_counts_mmr <- do.call(cbind, pb_list_mmr)

# 2f. Sample metadata table for DESeq2
pb_meta_mmr <- data.frame(
  sample    = names(pb_list_mmr),
  PID       = sub("__.*$", "", names(pb_list_mmr)),
  MMRStatus = sub("^.*__", "", names(pb_list_mmr)),
  row.names = names(pb_list_mmr),
  stringsAsFactors = FALSE
)
pb_meta_mmr$MMRStatus <- factor(pb_meta_mmr$MMRStatus, levels = c("MMRp", "MMRd"))
# MMRp is the reference level (denominator); MMRd is the numerator

cat(sprintf("  Pseudobulk samples: %d total (MMRd=%d, MMRp=%d)\n",
            nrow(pb_meta_mmr),
            sum(pb_meta_mmr$MMRStatus == "MMRd"),
            sum(pb_meta_mmr$MMRStatus == "MMRp")))

stopifnot(sum(pb_meta_mmr$MMRStatus == "MMRd") >= 3,
          sum(pb_meta_mmr$MMRStatus == "MMRp") >= 3)

# 2g. Filter lowly expressed genes
keep_genes_mmr <- rowSums(pb_counts_mmr >= 1) >= 2
pb_counts_mmr  <- pb_counts_mmr[keep_genes_mmr, ]

# 2h. DESeq2
dds_mmr <- DESeqDataSetFromMatrix(
  countData = as.matrix(pb_counts_mmr),
  colData   = pb_meta_mmr,
  design    = ~ MMRStatus
)
dds_mmr <- DESeq(dds_mmr, quiet = TRUE)

# 2i. Results: MMRd vs MMRp
de_mmr      <- as.data.frame(results(dds_mmr, contrast = c("MMRStatus", "MMRd", "MMRp")))
de_mmr$gene <- rownames(de_mmr)

write.csv(de_mmr, file.path(result_dir, "de_mmrd_vs_mmrp_myeloid.csv"), row.names = FALSE)
cat(sprintf("  DE genes (|log2FC| > 0.5, padj < 0.05): up in MMRd=%d, up in MMRp=%d\n",
            sum(de_mmr$log2FoldChange >  0.5 & de_mmr$padj < 0.05, na.rm = TRUE),
            sum(de_mmr$log2FoldChange < -0.5 & de_mmr$padj < 0.05, na.rm = TRUE)))

# 2j. Volcano plot
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

rm(seu_mmr, counts_mmr, meta_mmr, pb_list_mmr, pb_counts_mmr, pb_meta_mmr, dds_mmr); gc()


# =============================================================================
# STEP 3: Tumor vs Normal — per myeloid subcluster (pseudobulk DESeq2)
# Input:  seu,  cl295v11SubShort (cM01–cM10),  SPECIMEN_TYPE
#         Skips clusters with < 3 pseudobulk samples in either condition.
# Output: de_TvsN_<cluster>.csv  (one file per subcluster)
#         de_TvsN_all_subclusters.csv  (all clusters combined)
# =============================================================================
cat("\n--- STEP 3: Tumor vs Normal per myeloid subcluster ---\n")

subclusters    <- unique(seu$cl295v11SubShort)
de_per_cluster <- list()

for (sc in subclusters) {
  cat(sprintf("  Processing %s...\n", sc))

  # 3a. Subset to this subcluster's T and N cells
  seu_sc <- subset(seu, subset = cl295v11SubShort == sc &
                                 SPECIMEN_TYPE %in% c("T", "N"))

  # 3b. Raw counts for this subcluster
  counts_sc <- GetAssayData(seu_sc, slot = "counts", assay = "RNA")

  # 3c. Pseudobulk group labels
  meta_sc <- seu_sc@meta.data[, c("PID", "SPECIMEN_TYPE")]
  meta_sc$pb_group <- paste(meta_sc$PID, meta_sc$SPECIMEN_TYPE, sep = "__")

  # 3d. Sum counts per PID × condition group
  groups_sc <- unique(meta_sc$pb_group)
  pb_list_sc <- lapply(groups_sc, function(g) {
    cells_in_group <- rownames(meta_sc)[meta_sc$pb_group == g]
    if (length(cells_in_group) < 10) return(NULL)
    Matrix::rowSums(counts_sc[, cells_in_group, drop = FALSE])
  })
  names(pb_list_sc) <- groups_sc
  pb_list_sc <- pb_list_sc[!sapply(pb_list_sc, is.null)]

  pb_counts_sc <- do.call(cbind, pb_list_sc)

  # 3e. Sample metadata
  pb_meta_sc <- data.frame(
    PID           = sub("__.*$", "", names(pb_list_sc)),
    SPECIMEN_TYPE = sub("^.*__", "", names(pb_list_sc)),
    row.names     = names(pb_list_sc),
    stringsAsFactors = FALSE
  )
  pb_meta_sc$SPECIMEN_TYPE <- factor(pb_meta_sc$SPECIMEN_TYPE, levels = c("N", "T"))

  n_t <- sum(pb_meta_sc$SPECIMEN_TYPE == "T")
  n_n <- sum(pb_meta_sc$SPECIMEN_TYPE == "N")

  # 3f. Need at least 3 patients per condition to run DESeq2
  if (n_t < 3 | n_n < 3) {
    cat(sprintf("    Skipped (pseudobulk T=%d, N=%d patients)\n", n_t, n_n))
    next
  }

  # 3g. Filter genes, run DESeq2, extract results
  keep_sc      <- rowSums(pb_counts_sc >= 1) >= 2
  pb_counts_sc <- pb_counts_sc[keep_sc, ]

  de_sc <- tryCatch({
    dds_sc <- DESeqDataSetFromMatrix(
      countData = as.matrix(pb_counts_sc),
      colData   = pb_meta_sc,
      design    = ~ SPECIMEN_TYPE
    )
    dds_sc <- DESeq(dds_sc, quiet = TRUE)
    res    <- as.data.frame(results(dds_sc, contrast = c("SPECIMEN_TYPE", "T", "N")))
    res$gene    <- rownames(res)
    res$cluster <- sc
    res
  }, error = function(e) { message("    Error: ", e$message); NULL })

  if (!is.null(de_sc)) {
    de_per_cluster[[sc]] <- de_sc
    write.csv(de_sc, file.path(result_dir, paste0("de_TvsN_", sc, ".csv")),
              row.names = FALSE)
    cat(sprintf("    %d DE genes found (padj < 0.05)\n",
                sum(de_sc$padj < 0.05, na.rm = TRUE)))
  }
}

# 3h. Combine all subclusters into one table
if (length(de_per_cluster) > 0) {
  de_all <- bind_rows(de_per_cluster)
  write.csv(de_all, file.path(result_dir, "de_TvsN_all_subclusters.csv"),
            row.names = FALSE)
}


# =============================================================================
# STEP 4: Macrophage (cM02) vs Monocyte (cM01) — cell-type marker contrast
# Input:  seu subset to cM01 + cM02
# Output: de_macrophage_vs_monocyte.csv
#         volcano_macrophage_vs_monocyte.pdf / .png
# Note:   Wilcoxon is correct here — comparing two cell identities within the
#         same samples, not a condition across patients (no pseudoreplication).
# =============================================================================
cat("\n--- STEP 4: Macrophage (cM02) vs Monocyte (cM01) markers ---\n")

# 4a. Subset to just these two clusters
seu_mac_mono <- subset(seu, subset = cl295v11SubShort %in% c("cM02", "cM01"))

# 4b. Set cluster as the active identity so FindMarkers knows what to compare
Idents(seu_mac_mono) <- "cl295v11SubShort"

# 4c. Wilcoxon rank-sum test: cM02 (macrophage) vs cM01 (monocyte)
#     avg_log2FC > 0 means higher in macrophage
de_mac_mono       <- FindMarkers(seu_mac_mono, ident.1 = "cM02", ident.2 = "cM01",
                                 test.use = "wilcox", logfc.threshold = 0.1,
                                 min.pct = 0.1, verbose = FALSE)
de_mac_mono$gene  <- rownames(de_mac_mono)
write.csv(de_mac_mono, file.path(result_dir, "de_macrophage_vs_monocyte.csv"))

# 4d. Volcano plot (note: x-axis is avg_log2FC, y-axis is p_val_adj from Wilcoxon)
p_volcano_mac <- EnhancedVolcano(de_mac_mono,
  lab = de_mac_mono$gene, x = "avg_log2FC", y = "p_val_adj",
  title    = "Macrophage (cM02) vs Monocyte (cM01)",
  subtitle = "Wilcoxon (cell-type marker comparison)",
  pCutoff  = 0.05, FCcutoff = 0.5,
  pointSize = 1.5, labSize = 3,
  drawConnectors = TRUE, widthConnectors = 0.3)

ggsave(file.path(result_dir, "volcano_macrophage_vs_monocyte.pdf"), p_volcano_mac, width = 10, height = 8)
ggsave(file.path(result_dir, "volcano_macrophage_vs_monocyte.png"), p_volcano_mac, width = 10, height = 8, dpi = 200)


# =============================================================================
# STEP 5: Top DE genes heatmap
# Input:  de_tn  (Step 1 result, still in memory)
#         seu    (full myeloid object, for average expression per subcluster)
# Output: de_heatmap_TvsN_top_genes.pdf
# =============================================================================
cat("\n--- STEP 5: Top DE genes heatmap ---\n")

# 5a. Select top 30 up- and down-regulated genes from Step 1
top_up    <- de_tn %>% filter(padj < 0.05) %>% slice_max(n = 30, order_by = log2FoldChange)
top_down  <- de_tn %>% filter(padj < 0.05) %>% slice_min(n = 30, order_by = log2FoldChange)
top_genes <- c(top_up$gene, top_down$gene)
top_genes <- top_genes[top_genes %in% rownames(seu)]   # keep only genes in the object

# 5b. Compute average expression of these genes per myeloid subcluster
avg_expr <- AverageExpression(seu, features = top_genes, group.by = "cl295v11SubFull")
avg_mat  <- as.matrix(avg_expr$RNA)

# 5c. Z-score each gene across subclusters (so heatmap shows relative expression)
avg_mat_scaled <- t(scale(t(avg_mat)))

# 5d. Draw heatmap
pdf(file.path(result_dir, "de_heatmap_TvsN_top_genes.pdf"), width = 10, height = 12)
ht <- Heatmap(avg_mat_scaled, name = "Z-score",
              col             = colorRamp2(c(-2, 0, 2), c("#3498DB", "white", "#E74C3C")),
              row_names_gp    = gpar(fontsize = 7),
              column_names_gp = gpar(fontsize = 8),
              column_title    = "Top DE Genes: Tumor vs Normal (Pseudobulk DESeq2)")
draw(ht)
dev.off()

cat("\n Phase 5 complete!\n")
