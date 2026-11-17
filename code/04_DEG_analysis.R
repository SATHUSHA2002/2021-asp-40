
# ============================================================
# TB Gene Expression Analysis
# Script 04: Differential Expression Analysis
# ============================================================

source("code/02_preprocessing.R")

library(limma)

# 1. Prepare sample metadata
selected_samples$Group <- relevel(
  factor(selected_samples$Group),
  ref = "Control"
)

selected_samples$Gender <- factor(
  selected_samples$Gender
)

selected_samples$Ethnicity <- factor(
  selected_samples$Ethnicity
)

# 2. Create adjusted design matrix
design <- model.matrix(
  ~ Group + Gender + Ethnicity,
  data = selected_samples
)

cat("\nDesign matrix columns:\n")
print(colnames(design))

# 3. Check model is estimable
if (qr(design)$rank < ncol(design)) {
  stop("Design matrix is not full rank. Check confounding.")
}

# 4. Fit linear model
fit <- lmFit(
  selected_expression,
  design
)

fit <- eBayes(fit)

# 5. Extract TB vs Control results
deg_results <- topTable(
  fit,
  coef = "GroupTB",
  number = Inf,
  adjust.method = "BH",
  sort.by = "P"
)

# 6. Add probe IDs
deg_results$Probe_ID <- rownames(deg_results)

# 7. Select significant DEGs
significant_deg <- subset(
  deg_results,
  adj.P.Val < 0.05 & abs(logFC) > 1
)

# 8. Display summary
cat("\n===== DEG ANALYSIS SUMMARY =====\n")

cat("Total tested probes:",
    nrow(deg_results), "\n")

cat("Significant probes:",
    nrow(significant_deg), "\n")

cat("Upregulated probes:",
    sum(significant_deg$logFC > 1), "\n")

cat("Downregulated probes:",
    sum(significant_deg$logFC < -1), "\n")

# 9. Save results
write.csv(
  deg_results,
  "results/04_all_DEG_results.csv",
  row.names = FALSE
)

write.csv(
  significant_deg,
  "results/04_significant_DEG_results.csv",
  row.names = FALSE
)

cat("\nDEG analysis completed successfully!\n")
