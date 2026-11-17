
# ============================================================
# TB Gene Expression Analysis
# Script 07: KEGG Pathway Enrichment
# ============================================================

source("code/06_GO_enrichment.R")

library(clusterProfiler)
library(enrichplot)
library(ggplot2)

# 1. Use mapped significant genes
kegg_genes <- intersect(
  sig_entrez,
  background_entrez
)

# 2. Run KEGG enrichment
kegg_results <- enrichKEGG(
  gene = kegg_genes,
  universe = background_entrez,
  organism = "hsa",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

# 3. Convert KEGG IDs to gene symbols
if (nrow(as.data.frame(kegg_results)) > 0) {
  kegg_results <- setReadable(
    kegg_results,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID"
  )
}

# 4. Save results
kegg_table <- as.data.frame(kegg_results)

write.csv(
  kegg_table,
  "results/07_KEGG_pathway_results.csv",
  row.names = FALSE
)

# 5. Display summary
cat("\n===== KEGG ENRICHMENT SUMMARY =====\n")

cat("Input Entrez IDs:",
    length(kegg_genes), "\n")

cat("Enriched KEGG pathways:",
    nrow(kegg_table), "\n")

if (nrow(kegg_table) > 0) {
  
  cat("\nTop 10 KEGG pathways:\n")
  
  print(
    head(
      kegg_table[, c("Description", "p.adjust")],
      10
    )
  )
  
  # 6. Generate KEGG dotplot
  kegg_plot <- dotplot(
    kegg_results,
    showCategory = 20,
    title = "KEGG Pathway Enrichment - TB"
  )
  
  print(kegg_plot)
  
  ggsave(
    "figures/07_KEGG_dotplot.png",
    plot = kegg_plot,
    width = 11,
    height = 8,
    dpi = 300
  )
}

cat("\nKEGG enrichment analysis completed!\n")
