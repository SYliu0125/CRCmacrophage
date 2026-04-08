# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Reanalysis of Pelka et al. 2021 (Nature Immunology) — "Spatially organized multicellular immune hubs in human colorectal cancer". GEO accession: GSE178341. Focus is on myeloid/macrophage biology, with two key parallel comparisons: **Normal → pMMR** and **Normal → dMMR** (these are independent branches, not a linear progression).

## Running the Pipeline

All scripts must be run from the `analysis/` directory:

```bash
cd analysis/

# Install packages first (one-time)
Rscript code/00_install_packages.R

# Run full pipeline
Rscript code/run_all.R

# Run individual phases
Rscript code/01_load_and_qc.R
Rscript code/05_differential_expression.R
# etc.
```

Requires ≥32 GB RAM. Scripts are numbered 00–10 and must run sequentially — each phase reads `.rds` files produced by the previous phase.

## Architecture

### Data Flow (Sequential Dependencies)

```
analysis/data/           → 01_load_and_qc.R    → results/01_qc/seu_full_qc.rds
seu_full_qc.rds          → 02_processing.R     → results/02_processing/seu_full_processed.rds
seu_full_processed.rds   → 03_composition.R    → results/03_composition/ (plots only)
seu_full_processed.rds   → 04_myeloid_analysis.R → results/04_myeloid/seu_myeloid.rds
seu_myeloid.rds          → 05_differential_expression.R → results/05_de/
seu_myeloid.rds          → 06_pathway_analysis.R        → results/06_pathway/
seu_myeloid.rds          → 07_trajectory_analysis.R     → results/07_trajectory/
seu_myeloid.rds          → 08_regulon_scenic.R          → results/08_scenic/
```

Phases 05–08 all independently load `results/04_myeloid/seu_myeloid.rds`.

### Key Seurat Object Metadata Columns

- `cl295v11SubShort` — myeloid subcluster IDs (cM01–cM10); cM01 = monocytes (trajectory root), cM02 = macrophages
- `cl295v11SubFull` — full subcluster labels
- `clTopLevel` — 7 major cell types
- `SPECIMEN_TYPE` — `"T"` (tumor) or `"N"` (normal)
- `MMRStatus` — `"MMRp"` (proficient) or `"MMRd"` (deficient)
- `PID` — patient ID (used for pseudobulk aggregation)
- `Condition` — derived column created in script 07: `"N"`, `"MMRp"`, or `"MMRd"`

### DE Method Convention (Script 05)

- **Condition comparisons** (T vs N, MMRd vs MMRp): pseudobulk DESeq2 — counts aggregated per `PID × condition`, requires ≥3 pseudobulk samples per group
- **Cell-type marker contrasts** (cM02 vs cM01): Wilcoxon via `FindMarkers` — appropriate because these compare cell identities within samples, not conditions across patients

### SCENIC Fallback (Script 08)

If cisTarget databases are unavailable, script 08 falls back to scoring 19 curated myeloid TF regulons via AUCell. Key TFs: SPI1/PU.1, CEBPB, MAFB, MAF, STAT1, IRF1, RELA, HIF1A, TREM2.

## Condition Color Scheme

Used consistently across all visualization scripts:
- Normal: `#3498DB` (blue)
- MMRp: `#2ECC71` (green)  
- MMRd: `#E74C3C` (red)
