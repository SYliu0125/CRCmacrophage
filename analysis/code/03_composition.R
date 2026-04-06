# =============================================================================
# 03_composition.R
# Cell type composition analysis: proportions, Tumor vs Normal, MMRd vs MMRp
# =============================================================================

library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(cowplot)
library(viridis)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)

set.seed(42)
cat("=== Phase 3: Cell Type Composition Analysis ===\n\n")

result_dir <- "results/03_composition"
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

# --- 1. Load processed Seurat object ---
cat("Loading Seurat object...\n")
seu <- readRDS("results/02_processing/seu_full_processed.rds")
meta <- seu@meta.data

# --- 2. Cell type proportions per patient ---
cat("Computing cell type proportions per patient...\n")

prop_df <- meta %>%
  group_by(PID, clTopLevel) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(PID) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

# Add specimen info
patient_info <- meta %>%
  select(PID, SPECIMEN_TYPE, MMRStatus, TissueSiteSimple) %>%
  distinct()

prop_df <- left_join(prop_df, patient_info, by = "PID")

# --- 3. Stacked bar plot per patient ---
celltype_colors <- c(
  "Epi"    = "#E41A1C", "TNKILC" = "#377EB8", "Myeloid" = "#4DAF4A",
  "Plasma" = "#984EA3", "B"      = "#FF7F00", "Strom"   = "#A65628",
  "Mast"   = "#F781BF"
)

# Order patients by epithelial proportion
epi_order <- prop_df %>%
  filter(clTopLevel == "Epi") %>%
  arrange(desc(prop)) %>%
  pull(PID)

prop_df$PID <- factor(prop_df$PID, levels = epi_order)

p1 <- ggplot(prop_df, aes(x = PID, y = prop, fill = clTopLevel)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = celltype_colors) +
  labs(title = "Cell Type Proportions per Patient", x = "Patient", y = "Proportion", fill = "Cell Type") +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6))

ggsave(file.path(result_dir, "stacked_bar_by_patient.pdf"), p1, width = 14, height = 5)
ggsave(file.path(result_dir, "stacked_bar_by_patient.png"), p1, width = 14, height = 5, dpi = 150)

# --- 4. Tumor vs Normal comparison ---
cat("Comparing Tumor vs Normal...\n")

tn_prop <- meta %>%
  filter(SPECIMEN_TYPE %in% c("T", "N")) %>%
  group_by(SPECIMEN_TYPE, clTopLevel) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(SPECIMEN_TYPE) %>%
  mutate(prop = n / sum(n))

p2 <- ggplot(tn_prop, aes(x = clTopLevel, y = prop, fill = SPECIMEN_TYPE)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("T" = "#E74C3C", "N" = "#3498DB"),
                    labels = c("T" = "Tumor", "N" = "Normal")) +
  labs(title = "Cell Type Proportions: Tumor vs Normal",
       x = "Cell Type", y = "Proportion", fill = "Tissue") +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave(file.path(result_dir, "tumor_vs_normal.pdf"), p2, width = 8, height = 5)
ggsave(file.path(result_dir, "tumor_vs_normal.png"), p2, width = 8, height = 5, dpi = 150)

# Per-patient Tumor vs Normal (paired where available)
tn_patient <- meta %>%
  filter(SPECIMEN_TYPE %in% c("T", "N")) %>%
  group_by(PID, SPECIMEN_TYPE, clTopLevel) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(PID, SPECIMEN_TYPE) %>%
  mutate(prop = n / sum(n))

# Statistical test (Wilcoxon) for each cell type
cat("\n--- Tumor vs Normal Statistical Tests (Wilcoxon) ---\n")
for (ct in unique(tn_patient$clTopLevel)) {
  tumor_props <- tn_patient %>% filter(clTopLevel == ct, SPECIMEN_TYPE == "T") %>% pull(prop)
  normal_props <- tn_patient %>% filter(clTopLevel == ct, SPECIMEN_TYPE == "N") %>% pull(prop)
  if (length(tumor_props) > 2 & length(normal_props) > 2) {
    wt <- wilcox.test(tumor_props, normal_props)
    cat(sprintf("  %-10s p=%.4f (T: %.3f vs N: %.3f)\n", ct, wt$p.value,
                median(tumor_props), median(normal_props)))
  }
}

# --- 5. MMRd vs MMRp comparison ---
cat("\nComparing MMRd vs MMRp...\n")

mmr_prop <- meta %>%
  filter(MMRStatus %in% c("MMRd", "MMRp")) %>%
  group_by(MMRStatus, clTopLevel) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(MMRStatus) %>%
  mutate(prop = n / sum(n))

p3 <- ggplot(mmr_prop, aes(x = clTopLevel, y = prop, fill = MMRStatus)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("MMRd" = "#E74C3C", "MMRp" = "#2ECC71")) +
  labs(title = "Cell Type Proportions: MMRd vs MMRp",
       x = "Cell Type", y = "Proportion", fill = "MMR Status") +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave(file.path(result_dir, "mmrd_vs_mmrp.pdf"), p3, width = 8, height = 5)
ggsave(file.path(result_dir, "mmrd_vs_mmrp.png"), p3, width = 8, height = 5, dpi = 150)

# --- 6. Composition heatmap ---
cat("Generating composition heatmap...\n")

# Create patient x celltype proportion matrix
prop_mat <- prop_df %>%
  select(PID, clTopLevel, prop) %>%
  pivot_wider(names_from = clTopLevel, values_from = prop, values_fill = 0) %>%
  tibble::column_to_rownames("PID") %>%
  as.matrix()

# Annotation for patients
patient_ann <- patient_info %>%
  filter(PID %in% rownames(prop_mat)) %>%
  tibble::column_to_rownames("PID")

ha <- HeatmapAnnotation(
  SpecimenType = patient_ann[rownames(prop_mat), "SPECIMEN_TYPE"],
  MMRStatus    = patient_ann[rownames(prop_mat), "MMRStatus"],
  TissueSite   = patient_ann[rownames(prop_mat), "TissueSiteSimple"],
  col = list(
    SpecimenType = c("T" = "#E74C3C", "N" = "#3498DB"),
    MMRStatus    = c("MMRd" = "#E74C3C", "MMRp" = "#2ECC71"),
    TissueSite   = c("left" = "#8E44AD", "right" = "#F39C12", "rectum" = "#1ABC9C")
  ),
  annotation_name_side = "left"
)

pdf(file.path(result_dir, "composition_heatmap.pdf"), width = 10, height = 8)
ht <- Heatmap(t(prop_mat), name = "Proportion",
              col = colorRamp2(c(0, 0.25, 0.5), c("white", "orange", "red")),
              top_annotation = ha,
              cluster_columns = TRUE, cluster_rows = FALSE,
              column_title = "Cell Type Composition by Patient",
              show_column_names = FALSE,
              row_names_side = "left")
draw(ht)
dev.off()

# --- 7. Myeloid-focused composition ---
cat("\nMyeloid subcluster composition by condition...\n")

myeloid_comp <- meta %>%
  filter(clTopLevel == "Myeloid", SPECIMEN_TYPE %in% c("T", "N")) %>%
  group_by(SPECIMEN_TYPE, cl295v11SubFull) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(SPECIMEN_TYPE) %>%
  mutate(prop = n / sum(n))

p4 <- ggplot(myeloid_comp, aes(x = cl295v11SubFull, y = prop, fill = SPECIMEN_TYPE)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("T" = "#E74C3C", "N" = "#3498DB"),
                    labels = c("T" = "Tumor", "N" = "Normal")) +
  labs(title = "Myeloid Subcluster Proportions: Tumor vs Normal",
       x = "", y = "Proportion", fill = "Tissue") +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

ggsave(file.path(result_dir, "myeloid_composition_TvsN.pdf"), p4, width = 10, height = 5)
ggsave(file.path(result_dir, "myeloid_composition_TvsN.png"), p4, width = 10, height = 5, dpi = 150)

# Save composition tables
write.csv(prop_df, file.path(result_dir, "celltype_proportions_by_patient.csv"), row.names = FALSE)
write.csv(myeloid_comp, file.path(result_dir, "myeloid_composition_by_specimen.csv"), row.names = FALSE)

cat("\n Phase 3 complete!\n")
