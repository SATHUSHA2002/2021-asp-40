
# ============================================================
# TB Gene Expression Analysis
# Script 03: Quality Control Visualization
# ============================================================

# Load preprocessing script
source("code/02_preprocessing.R")

library(ggplot2)

# 1. Boxplot of gene expression values
png(
  "figures/03_expression_boxplot.png",
  width = 1600,
  height = 800,
  res = 150
)

boxplot(
  selected_expression,
  outline = FALSE,
  las = 2,
  col = ifelse(
    selected_samples$Group == "TB",
    "tomato",
    "skyblue"
  ),
  main = "Gene Expression Distribution - GSE42830",
  xlab = "Samples",
  ylab = "Expression Values",
  cex.axis = 0.4
)

dev.off()

cat("\nBoxplot saved successfully!\n")

# 2. PCA Analysis
# Samples are rows, genes are columns
pca_input <- t(selected_expression)

# Remove probes with zero variance
variable_probes <- apply(
  pca_input,
  2,
  var
) > 0

pca_input <- pca_input[, variable_probes, drop = FALSE]

pca_result <- prcomp(
  pca_input,
  center = TRUE,
  scale. = FALSE
)

# 3. Calculate explained variance
variance_explained <- 100 *
  pca_result$sdev^2 /
  sum(pca_result$sdev^2)

# 4. Prepare PCA data
pca_data <- data.frame(
  PC1 = pca_result$x[, 1],
  PC2 = pca_result$x[, 2],
  Group = selected_samples$Group
)

# 5. Create PCA plot
pca_plot <- ggplot(
  pca_data,
  aes(x = PC1, y = PC2, color = Group)
) +
  geom_point(size = 3, alpha = 0.8) +
  labs(
    title = "PCA - TB vs Control (GSE42830)",
    x = paste0(
      "PC1 (",
      round(variance_explained[1], 2),
      "%)"
    ),
    y = paste0(
      "PC2 (",
      round(variance_explained[2], 2),
      "%)"
    )
  ) +
  theme_minimal()

print(pca_plot)

ggsave(
  filename = "figures/03_PCA_plot.png",
  plot = pca_plot,
  width = 8,
  height = 6,
  dpi = 300
)

cat("\nPCA plot saved successfully!\n")

cat("\nPC1 variance:",
    round(variance_explained[1], 2), "%")

cat("\nPC2 variance:",
    round(variance_explained[2], 2), "%\n")

cat("\nQuality control visualization completed!\n")
