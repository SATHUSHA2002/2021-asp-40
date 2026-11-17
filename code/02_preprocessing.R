
# ============================================================
# TB Gene Expression Analysis
# Script 02: Data Preprocessing
# ============================================================

# Load the previous script
source("code/01_load_data.R")

# 1. Convert expression data to numeric matrix
selected_expression <- as.matrix(selected_expression)
storage.mode(selected_expression) <- "numeric"

# 2. Check missing and infinite values
missing_values <- sum(is.na(selected_expression))
infinite_values <- sum(is.infinite(selected_expression))

cat("\nMissing values:", missing_values)
cat("\nInfinite values:", infinite_values, "\n")

# 3. Check duplicate probe IDs
duplicate_probes <- sum(
  duplicated(rownames(selected_expression))
)

cat("Duplicate probe IDs:", duplicate_probes, "\n")

# 4. Verify sample order
stopifnot(
  identical(
    colnames(selected_expression),
    as.character(selected_samples$geo_accession)
  )
)

# 5. Inspect expression distribution
cat("\nExpression summary:\n")
print(summary(as.vector(selected_expression)))

# 6. Check sample metadata
cat("\nGroup distribution:\n")
print(table(selected_samples$Group))

cat("\nGender distribution:\n")
print(table(selected_samples$Gender, useNA = "ifany"))

cat("\nEthnicity distribution:\n")
print(table(selected_samples$Ethnicity, useNA = "ifany"))

# 7. Save preprocessing summary
preprocessing_summary <- data.frame(
  Dataset = "GSE42830",
  Probes = nrow(selected_expression),
  Samples = ncol(selected_expression),
  TB = sum(selected_samples$Group == "TB"),
  Control = sum(selected_samples$Group == "Control"),
  Missing_Values = missing_values,
  Infinite_Values = infinite_values,
  Duplicate_Probe_IDs = duplicate_probes
)

write.csv(
  preprocessing_summary,
  "results/02_preprocessing_summary.csv",
  row.names = FALSE
)

cat("\nPreprocessing checks completed!\n")
