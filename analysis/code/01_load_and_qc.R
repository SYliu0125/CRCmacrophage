# =============================================================================
# 01_load_and_qc.R
# Load full CRC scRNA-seq data (370K cells), merge metadata, QC
# Paper: Pelka et al. 2021
# =============================================================================

library(Seurat)
library(hdf5r)
library(ggplot2)
library(dplyr)
library(patchwork)
library(cowplot)
library(viridis)
library(data.table)

set.seed(42)
cat("=== Phase 1: Data Loading & QC ===\n\n")

# --- Paths ---
data_dir   <- "data"
result_dir <- "results/01_qc"
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

# --- 1. Load expression matrix ---
cat("Loading H5 expression matrix (this may take several minutes)...\n")
h5_file <- file.path(data_dir, "GSE178341_crc10x_full_c295v4_submit.h5")
counts <- Read10X_h5(h5_file)
cat(sprintf("  Expression matrix: %d genes x %d cells\n", nrow(counts), ncol(counts)))

# --- 2. Load annotations ---
cat("Loading cluster annotations...\n")
clusters <- fread(file.path(data_dir, "GSE178341_crc10x_full_c295v4_submit_cluster.csv"))

cat("Loading patient/sample metadata...\n")
metadata <- fread(file.path(data_dir, "GSE178341_crc10x_full_c295v4_submit_metatables.csv"))

# --- 3. Merge metadata ---
cat("Merging annotations...\n")
clusters_df <- as.data.frame(clusters)
rownames(clusters_df) <- clusters_df[[1]]
clusters_df[[1]] <- NULL

metadata_df <- as.data.frame(metadata)
rownames(metadata_df) <- metadata_df[[1]]
metadata_df[[1]] <- NULL

combined_meta <- merge(clusters_df, metadata_df, by = "row.names", all.x = TRUE)
rownames(combined_meta) <- combined_meta$Row.names
combined_meta$Row.names <- NULL

# Align with expression matrix
overlap <- intersect(colnames(counts), rownames(combined_meta))
cat(sprintf("  Overlapping cells: %d / %d\n", length(overlap), ncol(counts)))
counts <- counts[, overlap]
combined_meta <- combined_meta[overlap, ]

# --- 4. Create Seurat object ---
cat("Creating Seurat object...\n")
seu <- CreateSeuratObject(counts = counts, meta.data = combined_meta,
                          project = "CRC_Pelka2021", min.cells = 3, min.features = 200)
cat(sprintf("  Seurat object: %d genes x %d cells\n", nrow(seu), ncol(seu)))

rm(counts, combined_meta, clusters, metadata, clusters_df, metadata_df)
gc()

# --- 5. QC metrics ---
cat("Computing QC metrics...\n")
seu[["percent.mt"]]   <- PercentageFeatureSet(seu, pattern = "^MT-")
seu[["percent.ribo"]] <- PercentageFeatureSet(seu, pattern = "^RP[SL]")

cat("\n--- QC Summary ---\n")
cat(sprintf("  nFeature_RNA: median=%d, range=[%d, %d]\n",
            median(seu$nFeature_RNA), min(seu$nFeature_RNA), max(seu$nFeature_RNA)))
cat(sprintf("  nCount_RNA:   median=%d, range=[%d, %d]\n",
            median(seu$nCount_RNA), min(seu$nCount_RNA), max(seu$nCount_RNA)))
cat(sprintf("  percent.mt:   median=%.2f%%, range=[%.2f%%, %.2f%%]\n",
            median(seu$percent.mt), min(seu$percent.mt), max(seu$percent.mt)))

# --- 6. Cell type composition overview ---
cat("\n--- Cell Type Composition (clTopLevel) ---\n")
ct_counts <- sort(table(seu$clTopLevel), decreasing = TRUE)
for (ct in names(ct_counts)) {
  cat(sprintf("  %-10s %6d (%.1f%%)\n", ct, ct_counts[ct],
              100 * ct_counts[ct] / sum(ct_counts)))
}

# --- 7. QC Plots ---
cat("\nGenerating QC plots...\n")

# Violin by cell type
p1 <- VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
               group.by = "clTopLevel", pt.size = 0, ncol = 3) &
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10))

ggsave(file.path(result_dir, "qc_violin_by_celltype.pdf"), p1, width = 14, height = 5)
ggsave(file.path(result_dir, "qc_violin_by_celltype.png"), p1, width = 14, height = 5, dpi = 150)

# Scatter by cell type
p2 <- ggplot(seu@meta.data, aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt)) +
  geom_point(size = 0.05, alpha = 0.2) +
  scale_color_viridis(option = "magma", limits = c(0, 30), oob = scales::squish) +
  scale_x_log10() + scale_y_log10() +
  facet_wrap(~clTopLevel, nrow = 2) +
  labs(title = "QC: Counts vs Features", color = "% MT") +
  theme_cowplot() +
  theme(strip.text = element_text(size = 10))

ggsave(file.path(result_dir, "qc_scatter_by_celltype.pdf"), p2, width = 14, height = 8)
ggsave(file.path(result_dir, "qc_scatter_by_celltype.png"), p2, width = 14, height = 8, dpi = 150)

# Cell type bar plot
ct_df <- data.frame(CellType = names(ct_counts), Count = as.numeric(ct_counts))
ct_df$CellType <- factor(ct_df$CellType, levels = ct_df$CellType)

p3 <- ggplot(ct_df, aes(x = CellType, y = Count, fill = CellType)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = scales::comma(Count)), vjust = -0.3, size = 3.5) +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Cell Type Composition (All Cells)", x = "", y = "Number of Cells") +
  theme_cowplot() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 30, hjust = 1, size = 11))

ggsave(file.path(result_dir, "celltype_barplot.pdf"), p3, width = 8, height = 5)
ggsave(file.path(result_dir, "celltype_barplot.png"), p3, width = 8, height = 5, dpi = 150)

# --- 8. Clinical overview ---
cat("\n--- Specimen Type ---\n")
print(table(seu$SPECIMEN_TYPE))
cat("\n--- MMR Status ---\n")
print(table(seu$MMRStatus))
cat("\n--- Tissue Site ---\n")
print(table(seu$TissueSiteSimple))

# --- 9. Save ---
cat("\nSaving Seurat object...\n")
saveRDS(seu, file.path(result_dir, "seu_full_qc.rds"))
cat(sprintf("  Saved: %s (%.1f GB)\n", file.path(result_dir, "seu_full_qc.rds"),
            file.size(file.path(result_dir, "seu_full_qc.rds")) / 1e9))

cat("\n Phase 1 complete!\n")
