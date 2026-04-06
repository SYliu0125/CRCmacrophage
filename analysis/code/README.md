# Reanalysis of CRC Single-Cell Data (Pelka et al. 2021)

## Background

**Paper**: "Spatially organized multicellular immune hubs in human colorectal cancer" (Pelka et al., 2021, Nature Immunology)

**Data** (GEO: GSE178341):
- `GSE178341_crc10x_full_c295v4_submit.h5` — Expression matrix (~1.2 GB, ~370K cells × 43K genes)
- `GSE178341_crc10x_full_c295v4_submit_cluster.csv` — Cell cluster annotations (10 myeloid subclusters, 7 major cell types)
- `GSE178341_crc10x_full_c295v4_submit_metatables.csv` — Patient/clinical metadata

**Scope**: Full analysis of all cell types with **deep-dive into myeloid/macrophages**, including **trajectory analysis** and **regulon analysis (SCENIC)**. Key biological question: **how do macrophage populations change from Normal → pMMR and Normal → dMMR** (two parallel comparisons, not a linear progression). Results will be compared to published findings.

**Execution**: Scripts designed to run on a high-RAM machine (≥32 GB recommended).

---

## Analysis Pipeline (10 Phases)

### Phase 0: Package Installation
**Script**: `00_install_packages.R`
- Install all required CRAN, Bioconductor, and GitHub packages
- Includes: Seurat, harmony, monocle3, SCENIC, slingshot, etc.

---

### Phase 1: Data Loading & QC
**Script**: `01_load_and_qc.R`
- Load full H5 expression matrix (370K cells)
- Merge cluster annotations + patient metadata
- QC metrics: nFeature, nCount, percent.mt, percent.ribo
- QC violin plots by cell type, scatter plots
- Save full Seurat object (`seu_full_qc.rds`)

---

### Phase 2: Normalization, Dimensionality Reduction & Clustering
**Script**: `02_processing.R`
- SCTransform normalization
- PCA + elbow plot
- Harmony batch correction (across patients/batches)
- UMAP embedding
- Verify published clusters on UMAP (clTopLevel, cl295v11SubShort, cl295v11SubFull)
- UMAP colored by cell type, patient, tissue, MMR status
- Save processed object (`seu_full_processed.rds`)

---

### Phase 3: Cell Type Composition Analysis
**Script**: `03_composition.R`
- Cell type proportions per patient (stacked bar plots)
- Tumor vs Normal composition comparison
- MMRd vs MMRp composition comparison
- Statistical testing (Fisher's exact, Wilcoxon)
- Heatmap of frequencies across patients

---

### Phase 4: Myeloid / Macrophage Deep-Dive
**Script**: `04_myeloid_analysis.R`
- Subset myeloid cells (~42K: cM01-cM10)
- Re-normalize, PCA, UMAP on myeloid subset
- Macrophage subtype characterization (M1/M2/SPP1+/C1Q+)
- `FindAllMarkers` for each subcluster
- Dot plots & violin plots for key macrophage markers
- Myeloid composition by tissue type, MMR status
- Save myeloid Seurat object (`seu_myeloid.rds`)

---

### Phase 5: Differential Expression
**Script**: `05_differential_expression.R`
- Tumor vs Normal DE within myeloid cells (Wilcoxon + MAST)
- MMRd vs MMRp DE within myeloid subclusters
- Macrophage subtype contrasts (e.g., M1-like vs M2-like)
- Volcano plots (`EnhancedVolcano`)
- Heatmaps of top DE genes (`ComplexHeatmap`)
- Export DE tables as CSV

---

### Phase 6: Pathway & Gene Set Enrichment
**Script**: `06_pathway_analysis.R`
- GSEA with msigdbr Hallmark + GO gene sets
- Pathway enrichment via clusterProfiler (ORA + GSEA)
- Compare pathways between macrophage subtypes
- Focus: immune signaling, macrophage polarization, phagocytosis, antigen presentation
- Dot plots, bar plots, enrichment maps

---

### Phase 7: Trajectory Analysis
**Script**: `07_trajectory_analysis.R`

**Part A — Global trajectory (all myeloid)**:
- **Slingshot**: pseudotime trajectory on myeloid UMAP (monocyte → macrophage)
- Root at monocytes (cM01), trace to macrophage subtypes
- Pseudotime-colored UMAP with trajectory curves

**Part B — Two parallel condition comparisons**:

> **IMPORTANT**: Comparisons are **Normal → pMMR** and **Normal → dMMR** (two independent
> branches from normal tissue, NOT a sequential Normal → pMMR → dMMR progression).

- **Separate Slingshot** runs for Normal+pMMR cells and Normal+dMMR cells
- **Pseudotime density plots**: Normal vs pMMR, Normal vs dMMR (with Wilcoxon p-values)
- **Pseudotime box/violin plots** per comparison
- **UMAP by condition** with trajectory overlay per comparison
- **Macrophage subtype composition** along pseudotime, split by condition per comparison
- **Gene expression dynamics**: loess-smoothed expression of 20 key genes (SPP1, C1QA, VEGFA, TNF, IL1B, TREM2, etc.) along pseudotime for each comparison
- **Condition-specific trajectory heatmaps**: binned gene expression along pseudotime per condition

**Part C — tradeSeq**:
- Fit GAM to identify trajectory-associated genes
- Heatmap of top dynamic genes ordered by pseudotime

**Part D — Monocle3** (if available):
- Independent trajectory inference with `graph_test`
- Monocle3 pseudotime density by condition

---

### Phase 8: Regulon Analysis (SCENIC)
**Script**: `08_regulon_scenic.R`

**Part A — Regulon scoring**:
- Full SCENIC pipeline (GENIE3 → RcisTarget → AUCell) **if cisTarget databases available**
- Fallback: curated myeloid TF regulon scoring via AUCell (19 regulons)
- Key TFs: SPI1/PU.1, CEBPB, CEBPA, MAFB, MAF, STAT1, IRF1, NFkB/RELA, STAT3, HIF1A, MYC, PPARG, NR1H3/LXR, IRF8, IRF4, IRF7, TREM2, EGR1, NFAT

**Part B — Cluster-level analysis**:
- Regulon activity heatmap by myeloid subcluster
- RSS (Regulon Specificity Score) for cluster-specific TFs
- Regulon-based UMAP visualization

**Part C — Condition-specific regulon analysis**:

> **IMPORTANT**: Two parallel comparisons: **Normal → pMMR** and **Normal → dMMR**.

- **Cluster × condition regulon heatmap**: shows TF activation across all subclusters in each condition
- **Macrophage-specific ΔAUC bar plots**: which regulons gain/lose activity going Normal→pMMR and Normal→dMMR (side-by-side)
- **Violin plots**: key TF regulon scores in macrophages (cM02) split by N/pMMR/dMMR
- **Statistical testing**: Kruskal-Wallis + pairwise Wilcoxon for every regulon across conditions
- **CSV output**: `regulon_macrophage_condition_comparison.csv` with `Delta_pMMR_vs_Normal` and `Delta_dMMR_vs_Normal` columns

---

### Phase 9: Comparison to Published Findings
**Script**: `09_comparison_published.R`
- Reproduce key published figures (UMAP, cell type proportions)
- Compare our cluster assignments vs published
- Validate myeloid subtype markers match published Table/Figure
- Correlation of DE results with published gene lists
- Summary of concordance and novel findings

---

### Phase 10: Summary Figures
**Script**: `10_summary_figures.R`
- Publication-quality multi-panel figures
- Combined overview of all major results

---

## Required Packages

| Category           | Packages |
|--------------------|----------|
| **Core**           | Seurat, SeuratObject, hdf5r, Matrix, data.table |
| **Visualization**  | ggplot2, patchwork, cowplot, viridis, RColorBrewer, gridExtra, ComplexHeatmap, EnhancedVolcano, pheatmap |
| **Bioconductor**   | SingleCellExperiment, scater, scran, MAST, clusterProfiler, org.Hs.eg.db, enrichplot, msigdbr |
| **Batch Correction** | harmony |
| **Trajectory**     | monocle3, slingshot, tradeSeq |
| **Regulon/SCENIC** | SCENIC, AUCell, RcisTarget, GENIE3 |
| **Utilities**      | dplyr, tidyverse, reshape2, igraph, uwot, Rtsne, scales |

## How to Run

```bash
cd analysis/
Rscript code/run_all.R       # Run entire pipeline
```

Or run individual phases:
```bash
Rscript code/00_install_packages.R   # Install packages first
Rscript code/01_load_and_qc.R       # Then each phase sequentially
Rscript code/02_processing.R
# ... etc
```

## Directory Structure
```
analysis/
├── code/
│   ├── 00_install_packages.R
│   ├── 01_load_and_qc.R
│   ├── 02_processing.R
│   ├── 03_composition.R
│   ├── 04_myeloid_analysis.R
│   ├── 05_differential_expression.R
│   ├── 06_pathway_analysis.R
│   ├── 07_trajectory_analysis.R
│   ├── 08_regulon_scenic.R
│   ├── 09_comparison_published.R
│   ├── 10_summary_figures.R
│   ├── run_all.R
│   └── README.md              ← This file
├── data/
│   ├── GSE178341_crc10x_full_c295v4_submit.h5
│   ├── GSE178341_crc10x_full_c295v4_submit_cluster.csv
│   └── GSE178341_crc10x_full_c295v4_submit_metatables.csv
└── results/
    ├── 01_qc/
    ├── 02_processing/
    ├── 03_composition/
    ├── 04_myeloid/
    ├── 05_de/
    ├── 06_pathway/
    ├── 07_trajectory/
    ├── 08_scenic/
    ├── 09_comparison/
    └── 10_summary/
```

## Verification Plan

- Check Seurat object dimensions match expected counts
- Verify UMAP separates major cell types
- Confirm known markers (CD68, CD14 in myeloid; CD3E in T cells)
- Cross-reference myeloid subtype markers with published Table
- Trajectory root should start at monocytes
- SCENIC regulons should include known myeloid TFs (SPI1, CEBPB, etc.)
