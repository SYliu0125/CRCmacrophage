#!/usr/bin/env Rscript
# =============================================================================
# run_all.R
# Master script to run the complete CRC scRNA-seq reanalysis pipeline
# Run on a machine with >= 32 GB RAM
#
# Usage: Rscript code/run_all.R
# Or run individual phases: Rscript code/01_load_and_qc.R
# =============================================================================

cat("================================================================\n")
cat("  CRC scRNA-seq Reanalysis Pipeline\n")
cat("  Paper: Pelka et al. 2021 (Nature Immunology)\n")
cat("  Focus: Myeloid / Macrophage Analysis\n")
cat("================================================================\n\n")

start_time <- Sys.time()

# Set working directory (assumes run from analysis/ directory)
# setwd("path/to/analysis")

scripts <- c(
  "code/00_install_packages.R",
  "code/01_load_and_qc.R",
  "code/02_processing.R",
  "code/03_composition.R",
  "code/04_myeloid_analysis.R",
  "code/05_differential_expression.R",
  "code/06_pathway_analysis.R",
  "code/07_trajectory_analysis.R",
  "code/08_regulon_scenic.R",
  "code/09_comparison_published.R",
  "code/10_summary_figures.R"
)

for (script in scripts) {
  cat(sprintf("\n\n>>> Running: %s <<<\n", script))
  cat(sprintf("    Time: %s\n", Sys.time()))
  cat(paste(rep("=", 60), collapse = ""), "\n")
  
  tryCatch({
    source(script, local = FALSE)
    cat(sprintf("\n    Completed: %s\n", script))
  }, error = function(e) {
    cat(sprintf("\n    ERROR in %s: %s\n", script, e$message))
    cat("    Continuing with next script...\n")
  })
}

end_time <- Sys.time()
elapsed <- difftime(end_time, start_time, units = "mins")

cat("\n\n================================================================\n")
cat(sprintf("  Pipeline complete! Total time: %.1f minutes\n", elapsed))
cat("  Results saved in: results/\n")
cat("================================================================\n")
