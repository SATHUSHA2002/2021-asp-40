
# ============================================================
# TB Gene Expression Analysis
# Script 06: GO Biological Process Enrichment
# ============================================================

source("code/05_gene_annotation.R")

library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(enrichplot)

# 1. Convert significant gene symbols to Entrez IDs
sig_mapping <- bitr(
  unique(unique_deg$Gene_Symbol),
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

sig_entrez <- unique(sig_mapping$ENTREZID)

# 2. Prepare all measured, annotated genes
background_annotation <- annotation_data[
  annotation_data$ID %in% rownames(selected_expression),
  c("ID", "Symbol"),
  drop = FALSE
]

background_symbols <- unique(
  trimws(as.character(background_annotation$Symbol))
)

background_symbols <- background_symbols[
  !is.na(background_symbols) &
    background_symbols != ""
]

# 3. Convert background symbols to Entrez IDs
background_mapping <- bitr(
  background_symbols,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

background_entrez <- unique(background_mapping$ENTREZID)

# 4. Run GO enrichment
go_results <- enrichGO(
  gene = intersect(sig_entrez, background_entrez),
  universe = background_entrez,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

# 5. Save results
go_table <- as.data.frame(go_results)

write.csv(
  go_table,
  "results/06_GO_BP_results.csv",
  row.names = FALSE
)

# 6. Display summary
cat("\n===== GO ENRICHMENT SUMMARY =====\n")

cat("Significant gene Entrez IDs:",
    length(sig_entrez), "\n")

cat("Background Entrez IDs:",
    length(background_entrez), "\n")

cat("Enriched GO BP terms:",
    nrow(go_table), "\n")

cat("\nTop GO Biological Processes:\n")
print(
  head(
    go_table[, c("Description", "p.adjust")],
    10
  )
)

# 7. Create GO dotplot
if (nrow(go_table) > 0) {
  
  go_plot <- dotplot(
    go_results,
    showCategory = 20,
    title = "GO Biological Process Enrichment"
  )
  
  print(go_plot)
  
  ggsave(
    "figures/06_GO_BP_dotplot.png",
    plot = go_plot,
    width = 11,
    height = 8,
    dpi = 300
  )
}

cat("\nGO enrichment analysis completed!\n")
