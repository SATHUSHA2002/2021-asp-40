
# ============================================================
# TB Gene Expression Analysis
# Script 09: Final Research Summary
# ============================================================

# Run the analysis scripts in sequence
source("code/03_quality_control.R")
source("code/07_KEGG_enrichment.R")
source("code/08_STRING_network.R")

# Create summary table
research_summary <- data.frame(
  Metric = c(
    "Dataset",
    "Total Samples",
    "TB Samples",
    "Control Samples",
    "Total Probes",
    "Significant Probes",
    "Upregulated Probes",
    "Downregulated Probes",
    "Annotated Significant Probes",
    "Unique Gene Symbols",
    "Mapped Significant Entrez IDs",
    "GO BP Terms",
    "KEGG Pathways",
    "STRING Proteins",
    "STRING Interactions"
  ),
  Value = c(
    "GSE42830",
    ncol(selected_expression),
    sum(selected_samples$Group == "TB"),
    sum(selected_samples$Group == "Control"),
    nrow(selected_expression),
    nrow(significant_deg),
    sum(significant_deg$logFC > 1),
    sum(significant_deg$logFC < -1),
    nrow(annotated_deg),
    nrow(unique_deg),
    length(sig_entrez),
    nrow(go_table),
    nrow(kegg_table),
    igraph::vcount(network),
    igraph::ecount(network)
  )
)

# Save final summary
write.csv(
  research_summary,
  "results/09_final_research_summary.csv",
  row.names = FALSE
)

cat("\n===== FINAL RESEARCH SUMMARY =====\n")
print(research_summary, row.names = FALSE)

cat("\nFinal research summary saved successfully!\n")
