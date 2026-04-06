# =============================================================================
# 06_pathway_analysis.R
# Pathway & Gene Set Enrichment Analysis on myeloid DE results
# =============================================================================

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(msigdbr)
library(ggplot2)
library(dplyr)
library(cowplot)
library(ComplexHeatmap)
library(circlize)

set.seed(42)
cat("=== Phase 6: Pathway & Gene Set Enrichment Analysis ===\n\n")

result_dir <- "results/06_pathway"
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

# --- 1. Load DE results ---
cat("Loading DE results...\n")
de_tn  <- read.csv("results/05_de/de_tumor_vs_normal_myeloid.csv", row.names = 1)
de_mmr <- read.csv("results/05_de/de_mmrd_vs_mmrp_myeloid.csv", row.names = 1)

# --- 2. Prepare gene lists ---
# For GSEA: ranked gene list by log2FC * -log10(p_val_adj)
prepare_ranked_list <- function(de_df) {
  de_df <- de_df %>%
    filter(!is.na(p_val_adj) & p_val_adj > 0) %>%
    mutate(rank_metric = avg_log2FC * -log10(p_val_adj)) %>%
    arrange(desc(rank_metric))
  
  ranked <- de_df$rank_metric
  names(ranked) <- de_df$gene
  return(ranked)
}

ranked_tn  <- prepare_ranked_list(de_tn)
ranked_mmr <- prepare_ranked_list(de_mmr)

# --- 3. Get MSigDB gene sets ---
cat("Loading MSigDB gene sets...\n")

# Hallmark
hallmark <- msigdbr(species = "Homo sapiens", category = "H") %>%
  select(gs_name, gene_symbol) %>%
  rename(term = gs_name, gene = gene_symbol)

# GO Biological Process
gobp <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP") %>%
  select(gs_name, gene_symbol) %>%
  rename(term = gs_name, gene = gene_symbol)

# Reactome
reactome <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME") %>%
  select(gs_name, gene_symbol) %>%
  rename(term = gs_name, gene = gene_symbol)

# =========================================================================
# GSEA: Tumor vs Normal
# =========================================================================
cat("\n--- GSEA: Tumor vs Normal ---\n")

# Hallmark GSEA
gsea_hallmark_tn <- GSEA(ranked_tn, TERM2GENE = hallmark, pvalueCutoff = 0.05,
                          minGSSize = 15, maxGSSize = 500, verbose = FALSE)

if (nrow(gsea_hallmark_tn@result) > 0) {
  write.csv(gsea_hallmark_tn@result, file.path(result_dir, "gsea_hallmark_TvsN.csv"), row.names = FALSE)
  
  p1 <- dotplot(gsea_hallmark_tn, showCategory = 20, split = ".sign") +
    facet_wrap(~.sign, scales = "free_y") +
    labs(title = "Hallmark GSEA: Tumor vs Normal (Myeloid)") +
    theme_cowplot() +
    theme(axis.text.y = element_text(size = 7))
  
  ggsave(file.path(result_dir, "gsea_hallmark_TvsN_dotplot.pdf"), p1, width = 12, height = 8)
  ggsave(file.path(result_dir, "gsea_hallmark_TvsN_dotplot.png"), p1, width = 12, height = 8, dpi = 200)
  
  # GSEA enrichment plot for top pathways
  top_paths <- head(gsea_hallmark_tn@result$ID, 4)
  pdf(file.path(result_dir, "gsea_hallmark_TvsN_enrichplots.pdf"), width = 10, height = 10)
  for (path in top_paths) {
    print(gseaplot2(gsea_hallmark_tn, geneSetID = path, title = path))
  }
  dev.off()
}

# GO BP GSEA
gsea_gobp_tn <- GSEA(ranked_tn, TERM2GENE = gobp, pvalueCutoff = 0.05,
                      minGSSize = 15, maxGSSize = 500, verbose = FALSE)

if (nrow(gsea_gobp_tn@result) > 0) {
  write.csv(gsea_gobp_tn@result, file.path(result_dir, "gsea_gobp_TvsN.csv"), row.names = FALSE)
  
  p2 <- dotplot(gsea_gobp_tn, showCategory = 20, split = ".sign") +
    facet_wrap(~.sign, scales = "free_y") +
    labs(title = "GO BP GSEA: Tumor vs Normal (Myeloid)") +
    theme_cowplot() +
    theme(axis.text.y = element_text(size = 6))
  
  ggsave(file.path(result_dir, "gsea_gobp_TvsN_dotplot.pdf"), p2, width = 14, height = 10)
  ggsave(file.path(result_dir, "gsea_gobp_TvsN_dotplot.png"), p2, width = 14, height = 10, dpi = 200)
}

# =========================================================================
# GSEA: MMRd vs MMRp
# =========================================================================
cat("\n--- GSEA: MMRd vs MMRp ---\n")

gsea_hallmark_mmr <- GSEA(ranked_mmr, TERM2GENE = hallmark, pvalueCutoff = 0.05,
                           minGSSize = 15, maxGSSize = 500, verbose = FALSE)

if (nrow(gsea_hallmark_mmr@result) > 0) {
  write.csv(gsea_hallmark_mmr@result, file.path(result_dir, "gsea_hallmark_MMR.csv"), row.names = FALSE)
  
  p3 <- dotplot(gsea_hallmark_mmr, showCategory = 20, split = ".sign") +
    facet_wrap(~.sign, scales = "free_y") +
    labs(title = "Hallmark GSEA: MMRd vs MMRp (Myeloid)") +
    theme_cowplot() +
    theme(axis.text.y = element_text(size = 7))
  
  ggsave(file.path(result_dir, "gsea_hallmark_MMR_dotplot.pdf"), p3, width = 12, height = 8)
  ggsave(file.path(result_dir, "gsea_hallmark_MMR_dotplot.png"), p3, width = 12, height = 8, dpi = 200)
}

# =========================================================================
# ORA: Over-representation analysis on significant DE genes
# =========================================================================
cat("\n--- ORA: Significant genes from Tumor vs Normal ---\n")

sig_up   <- de_tn %>% filter(p_val_adj < 0.05, avg_log2FC > 0.5) %>% pull(gene)
sig_down <- de_tn %>% filter(p_val_adj < 0.05, avg_log2FC < -0.5) %>% pull(gene)
background <- de_tn$gene

cat(sprintf("  Upregulated: %d genes, Downregulated: %d genes\n",
            length(sig_up), length(sig_down)))

# GO enrichment for upregulated genes
if (length(sig_up) > 10) {
  ego_up <- enrichGO(gene = sig_up, universe = background,
                     OrgDb = org.Hs.eg.db, keyType = "SYMBOL",
                     ont = "BP", pAdjustMethod = "BH",
                     pvalueCutoff = 0.05, readable = FALSE)
  
  if (nrow(ego_up@result) > 0) {
    write.csv(ego_up@result, file.path(result_dir, "ora_go_upregulated_TvsN.csv"), row.names = FALSE)
    
    p4 <- dotplot(ego_up, showCategory = 20) +
      labs(title = "GO BP: Upregulated in Tumor (Myeloid)") +
      theme_cowplot() +
      theme(axis.text.y = element_text(size = 7))
    
    ggsave(file.path(result_dir, "ora_go_up_TvsN.pdf"), p4, width = 10, height = 8)
    ggsave(file.path(result_dir, "ora_go_up_TvsN.png"), p4, width = 10, height = 8, dpi = 200)
  }
}

# GO enrichment for downregulated genes
if (length(sig_down) > 10) {
  ego_down <- enrichGO(gene = sig_down, universe = background,
                       OrgDb = org.Hs.eg.db, keyType = "SYMBOL",
                       ont = "BP", pAdjustMethod = "BH",
                       pvalueCutoff = 0.05, readable = FALSE)
  
  if (nrow(ego_down@result) > 0) {
    write.csv(ego_down@result, file.path(result_dir, "ora_go_downregulated_TvsN.csv"), row.names = FALSE)
    
    p5 <- dotplot(ego_down, showCategory = 20) +
      labs(title = "GO BP: Downregulated in Tumor (Myeloid)") +
      theme_cowplot() +
      theme(axis.text.y = element_text(size = 7))
    
    ggsave(file.path(result_dir, "ora_go_down_TvsN.pdf"), p5, width = 10, height = 8)
    ggsave(file.path(result_dir, "ora_go_down_TvsN.png"), p5, width = 10, height = 8, dpi = 200)
  }
}

cat("\n Phase 6 complete!\n")
