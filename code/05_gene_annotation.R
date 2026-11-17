
# ============================================================
# TB Gene Expression Analysis
# Script 05: Gene Annotation
# ============================================================

source("code/04_DEG_analysis.R")

library(GEOquery)

# 1. Download GPL10558 platform annotation
gpl <- getGEO(
  "GPL10558",
  destdir = "data"
)

annotation_data <- Table(gpl)

cat("\nAnnotation columns:\n")
print(colnames(annotation_data))

# 2. Inspect platform annotation
cat("\nFirst annotation rows:\n")
print(head(annotation_data))

# 3. Check required annotation fields
required_columns <- c(
  "ID",
  "Symbol"
)

if (!all(required_columns %in% colnames(annotation_data))) {
  stop(
    "Expected annotation columns not found. ",
    "Inspect colnames(annotation_data) before continuing."
  )
}

# 4. Prepare gene annotation table
gene_annotation <- annotation_data[
  ,
  required_columns,
  drop = FALSE
]

colnames(gene_annotation) <- c(
  "Probe_ID",
  "Gene_Symbol"
)

# 5. Match annotation with significant probes
annotated_deg <- merge(
  significant_deg,
  gene_annotation,
  by = "Probe_ID",
  all.x = TRUE,
  sort = FALSE
)

# 6. Remove missing gene symbols
annotated_deg <- annotated_deg[
  !is.na(annotated_deg$Gene_Symbol) &
    trimws(annotated_deg$Gene_Symbol) != "",
  ,
  drop = FALSE
]

# 7. Generate unique gene list
unique_deg <- annotated_deg[
  order(annotated_deg$adj.P.Val),
  ,
  drop = FALSE
]

unique_deg <- unique_deg[
  !duplicated(unique_deg$Gene_Symbol),
  ,
  drop = FALSE
]

# 8. Display results
cat("\n===== GENE ANNOTATION SUMMARY =====\n")

cat("Significant probes:",
    nrow(significant_deg), "\n")

cat("Annotated probes:",
    nrow(annotated_deg), "\n")

cat("Unique gene symbols:",
    nrow(unique_deg), "\n")

cat("\nTop 10 genes:\n")

print(
  head(
    unique_deg[
      ,
      c("Gene_Symbol", "logFC", "adj.P.Val")
    ],
    10
  )
)

# 9. Save annotation results
write.csv(
  annotated_deg,
  "results/05_annotated_DEG_results.csv",
  row.names = FALSE
)

write.csv(
  unique_deg,
  "results/05_unique_DEG_genes.csv",
  row.names = FALSE
)

cat("\nGene annotation completed!\n")
