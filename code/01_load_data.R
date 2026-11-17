
# ============================================================
# TB Gene Expression Analysis
# Script 01: Load GSE42830 Dataset
# ============================================================

# Load required package
library(GEOquery)
library(Biobase)

# Create output folders if necessary
dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

# 1. Load GEO Series Matrix from local data folder
data_file <- "data/GSE42830_series_matrix.txt.gz"

if (!file.exists(data_file)) {
  stop("GSE42830 dataset not found in data folder.")
}

gse <- getGEO(
  filename = data_file,
  GSEMatrix = TRUE
)

# 2. Extract expression matrix
expression_matrix <- exprs(gse)

# 3. Extract sample metadata
sample_metadata <- pData(gse)

# 4. Identify TB and healthy control samples
disease_labels <- trimws(
  as.character(sample_metadata[["disease state:ch1"]])
)

selected_samples <- sample_metadata[
  !is.na(disease_labels) &
    disease_labels %in% c("TB", "Control"),
  ,
  drop = FALSE
]

# 5. Create group labels
selected_samples$Group <- factor(
  selected_samples[["disease state:ch1"]],
  levels = c("Control", "TB")
)

# 6. Extract gender and ethnicity
selected_samples$Gender <- factor(
  trimws(as.character(
    selected_samples[["gender:ch1"]]
  ))
)

selected_samples$Ethnicity <- factor(
  trimws(as.character(
    selected_samples[["ethnicity:ch1"]]
  ))
)

# 7. Select expression data for TB and Control
selected_expression <- expression_matrix[
  ,
  as.character(selected_samples$geo_accession),
  drop = FALSE
]

# 8. Verify sample alignment
stopifnot(
  identical(
    colnames(selected_expression),
    as.character(selected_samples$geo_accession)
  )
)

# 9. Display dataset summary
cat("\n===== DATASET SUMMARY =====\n")
cat("Dataset: GSE42830\n")
cat("Total probes:", nrow(selected_expression), "\n")
cat("Selected samples:", ncol(selected_expression), "\n")

cat("\nSample groups:\n")
print(table(selected_samples$Group))

cat("\nExpression dimensions:\n")
print(dim(selected_expression))

cat("\nMissing expression values:",
    sum(is.na(selected_expression)), "\n")

cat("\nDataset loaded successfully!\n")
