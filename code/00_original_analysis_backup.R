# ============================================================
# TB Gene Expression Analysis Project
# Dataset: GSE42830
# Step 1: Load and inspect the GEO Series Matrix file
# ============================================================

# Check the files available in the data folder
list.files("data")

# Define the path to the dataset
data_file <- "data/GSE42830_series_matrix.txt.gz"

# Read the GEO Series Matrix file
# GEO matrix files contain metadata at the beginning.
# We first read all lines to find where the expression table starts.
all_lines <- readLines(gzfile(data_file))

# Find the line containing the expression data table
table_start <- grep("!series_matrix_table_begin", all_lines)

# Find the line where the expression data table ends
table_end <- grep("!series_matrix_table_end", all_lines)

# Read the expression data table
tb_expression <- read.delim(
  gzfile(data_file),
  skip = table_start,
  header = TRUE,
  sep = "\t",
  check.names = FALSE
)

# Remove the final GEO table end line if it was read
tb_expression <- tb_expression[
  !grepl("!series_matrix_table_end", tb_expression[, 1]),
]

# Display the first six rows of expression data
head(tb_expression)

# Check the number of probes and samples
dim(tb_expression)

# Display the first few sample IDs
colnames(tb_expression)[1:10]

# Step 2: Extract sample information from GEO metadata

# Find sample IDs, titles and characteristics
sample_info_lines <- all_lines[
  grepl(
    "^!Sample_(geo_accession|title|source_name_ch1|characteristics_ch1)",
    all_lines
  )
]

# Display sample information
cat(sample_info_lines, sep = "\n")


# Step 3: Create a sample information table

# Function to extract values from GEO metadata
get_values <- function(pattern) {
  line <- grep(pattern, all_lines, value = TRUE)
  
  values <- strsplit(line, "\t", fixed = TRUE)[[1]][-1]
  
  gsub('^"|"$', "", values)
}

# Extract sample IDs
sample_ids <- get_values("^!Sample_geo_accession\\t")

# Extract sample titles
sample_titles <- get_values("^!Sample_title\\t")

# Extract disease information
disease_line <- grep(
  '^!Sample_characteristics_ch1\\t"disease state:',
  all_lines,
  value = TRUE
)

disease_states <- strsplit(
  disease_line,
  "\t",
  fixed = TRUE
)[[1]][-1]

disease_states <- gsub('^"|"$', "", disease_states)
disease_states <- sub("disease state: ", "", disease_states)

# Combine information into a table
sample_info <- data.frame(
  Sample_ID = sample_ids,
  Sample_Title = sample_titles,
  Disease = disease_states
)

# Display the first 10 samples
head(sample_info, 10)

# Count the samples in each disease group
table(sample_info$Disease)

# Check whether sample information matches expression data
nrow(sample_info)
ncol(tb_expression) - 1

# Step 4: Select TB and healthy control samples

# Select only TB and Control samples
selected_samples <- sample_info[
  sample_info$Disease %in% c("TB", "Control"),
]

# Create a group column
selected_samples$Group <- ifelse(
  selected_samples$Disease == "TB",
  "TB",
  "Control"
)

# Set Control as the reference group
selected_samples$Group <- factor(
  selected_samples$Group,
  levels = c("Control", "TB")
)

# Select matching expression columns
selected_expression <- tb_expression[
  ,
  c("ID_REF", selected_samples$Sample_ID),
  drop = FALSE
]

# Check the selected sample information
head(selected_samples)

# Count samples in each group
table(selected_samples$Group)

# Check the dimensions of the selected expression data
dim(selected_expression)

# Verify sample IDs match
identical(
  colnames(selected_expression)[-1],
  selected_samples$Sample_ID
)

# Step 5: Inspect expression values

# Remove the probe ID column for numerical inspection
expression_values <- selected_expression[, -1]

# Check the data type of expression values
class(expression_values[[1]])

# Check the first few expression values
head(expression_values[, 1:5])

# Find the minimum and maximum values
range(
  as.matrix(expression_values),
  na.rm = TRUE
)

# Calculate summary statistics
summary(
  as.vector(as.matrix(expression_values))
)

# Count missing values
sum(is.na(expression_values))

# Count infinite values
sum(!is.finite(as.matrix(expression_values)))

# Step 6: Quality control using a boxplot

# Convert expression data into a numerical matrix
expression_matrix <- as.matrix(expression_values)

# Draw a boxplot for all 54 samples
boxplot(
  expression_matrix,
  outline = FALSE,
  las = 2,
  cex.axis = 0.6,
  main = "Expression Distribution Across Samples",
  xlab = "Samples",
  ylab = "Normalized Log2 Expression Values"
)

# Add a horizontal reference line at zero
abline(h = 0, col = "red", lty = 2)

# Step 7: Save the expression distribution boxplot

png(
  filename = "figures/01_expression_boxplot.png",
  width = 1400,
  height = 700,
  res = 120
)

boxplot(
  expression_matrix,
  outline = FALSE,
  las = 2,
  cex.axis = 0.6,
  main = "Expression Distribution Across Samples",
  xlab = "Samples",
  ylab = "Normalized Log2 Expression Values"
)

abline(h = 0, col = "red", lty = 2)

dev.off()

# Step 8: Principal Component Analysis (PCA)

# Transpose the expression matrix
# Rows = samples, columns = probes
pca_data <- t(expression_matrix)

# Perform PCA
pca_result <- prcomp(
  pca_data,
  center = TRUE,
  scale. = FALSE
)

# Calculate explained variance
variance_explained <- 100 * (
  pca_result$sdev^2 / sum(pca_result$sdev^2)
)

# Extract PCA coordinates
pca_df <- data.frame(
  Sample_ID = rownames(pca_result$x),
  PC1 = pca_result$x[, 1],
  PC2 = pca_result$x[, 2]
)

# Match sample groups to PCA results
pca_df$Group <- selected_samples$Group[
  match(pca_df$Sample_ID, selected_samples$Sample_ID)
]

# Assign colours to groups
group_colors <- c(
  "Control" = "blue",
  "TB" = "red"
)

# Draw PCA plot
plot(
  pca_df$PC1,
  pca_df$PC2,
  col = group_colors[as.character(pca_df$Group)],
  pch = 19,
  xlab = paste0("PC1 (", round(variance_explained[1], 2), "%)"),
  ylab = paste0("PC2 (", round(variance_explained[2], 2), "%)"),
  main = "PCA: TB vs Healthy Controls"
)

# Add legend
legend(
  "topright",
  legend = c("Control", "TB"),
  col = c("blue", "red"),
  pch = 19
)

# Step 9: Save PCA plot

png(
  filename = "figures/02_PCA_TB_vs_Control.png",
  width = 1200,
  height = 800,
  res = 120
)

plot(
  pca_df$PC1,
  pca_df$PC2,
  col = group_colors[as.character(pca_df$Group)],
  pch = 19,
  xlab = paste0("PC1 (", round(variance_explained[1], 2), "%)"),
  ylab = paste0("PC2 (", round(variance_explained[2], 2), "%)"),
  main = "PCA: TB vs Healthy Controls"
)

legend(
  "topright",
  legend = c("Control", "TB"),
  col = c("blue", "red"),
  pch = 19
)

dev.off()

# Step 10C: Create the design matrix

# Load limma
library(limma)

# Create design matrix using sample groups
design <- model.matrix(
  ~ 0 + Group,
  data = selected_samples
)

# Rename columns
colnames(design) <- c("Control", "TB")

# Display the first 10 rows
head(design, 10)

# Check design matrix dimensions
dim(design)

# Check number of samples in each group
colSums(design)

# Verify expression samples match the design
identical(
  colnames(expression_matrix),
  selected_samples$Sample_ID
)
# Step 10D: Fit the differential expression model

# Define the comparison: TB versus Control
contrast_matrix <- makeContrasts(
  TB_vs_Control = TB - Control,
  levels = design
)

# Display the contrast matrix
contrast_matrix

# Fit a linear model for each probe
fit <- lmFit(expression_matrix, design)

# Apply the TB versus Control contrast
fit2 <- contrasts.fit(fit, contrast_matrix)

# Apply empirical Bayes moderation
fit2 <- eBayes(fit2)

# Display a summary of differential expression
summary(decideTests(fit2))

# Step 10E: Inspect differential expression results

# Extract results for all probes
deg_results <- topTable(
  fit2,
  coef = "TB_vs_Control",
  number = Inf,
  adjust.method = "BH"
)

# Display first 10 results
head(deg_results, 10)

# Check the columns
colnames(deg_results)

# Summary of estimated expression differences
summary(deg_results$logFC)

# Summary of adjusted p-values
summary(deg_results$adj.P.Val)

# Count probes with large differences and FDR < 0.05
sum(
  abs(deg_results$logFC) > 1 &
    deg_results$adj.P.Val < 0.05,
  na.rm = TRUE
)

# Step 10F: Verify probe identifiers

# Check original probe IDs
head(selected_expression$ID_REF, 10)

# Check whether probe IDs are unique
anyDuplicated(selected_expression$ID_REF)

# Check the result row names
head(rownames(deg_results), 10)

# Check number of results
nrow(deg_results)

# Check whether expression matrix has probe IDs
head(rownames(expression_matrix), 10)

# Step 10G: Correct probe identifiers

# Assign original probe IDs to the expression matrix
rownames(expression_matrix) <- selected_expression$ID_REF

# Verify the corrected probe IDs
head(rownames(expression_matrix), 10)

# Repeat limma analysis using corrected matrix
fit <- lmFit(expression_matrix, design)

fit2 <- contrasts.fit(fit, contrast_matrix)

fit2 <- eBayes(fit2)

# Extract corrected differential expression results
deg_results <- topTable(
  fit2,
  coef = "TB_vs_Control",
  number = Inf,
  adjust.method = "BH"
)

# Display the first 10 results
head(deg_results, 10)

# Verify total number of probes
nrow(deg_results)

# Verify significance counts
summary(decideTests(fit2))

# Step 11: Save differential expression results

# Convert row names into a Probe_ID column
deg_results$Probe_ID <- rownames(deg_results)

# Move Probe_ID to the first column
deg_results <- deg_results[
  ,
  c("Probe_ID", "logFC", "AveExpr",
    "t", "P.Value", "adj.P.Val", "B")
]

# Save all probe results
write.csv(
  deg_results,
  file = "results/01_all_DEG_results.csv",
  row.names = FALSE
)

# Select significant probes with larger differences
significant_deg <- deg_results[
  deg_results$adj.P.Val < 0.05 &
    abs(deg_results$logFC) > 1,
]

# Save selected significant probes
write.csv(
  significant_deg,
  file = "results/02_significant_DEG_results.csv",
  row.names = FALSE
)

# Check number of selected probes
nrow(significant_deg)

# Check saved files
list.files("results")

# Step 12: Compare TB and Control expression distributions

# Calculate the median expression for each sample
sample_medians <- apply(
  expression_matrix,
  2,
  median
)

# Create a table of sample medians and groups
median_data <- data.frame(
  Sample_ID = colnames(expression_matrix),
  Group = selected_samples$Group,
  Median_Expression = sample_medians
)

# Summarize medians by group
aggregate(
  Median_Expression ~ Group,
  data = median_data,
  FUN = summary
)

# Draw a boxplot comparing the groups
boxplot(
  Median_Expression ~ Group,
  data = median_data,
  col = c("lightblue", "salmon"),
  main = "Sample Median Expression: TB vs Control",
  xlab = "Group",
  ylab = "Median Expression"
)

# Find GEO data processing information
processing_lines <- grep(
  "^!Sample_data_processing",
  all_lines,
  value = TRUE
)

# Display the first processing entry
if (length(processing_lines) > 0) {
  cat(
    strsplit(processing_lines[1], "\t")[[1]][1:2],
    sep = "\n"
  )
} else {
  cat("No sample processing information found.\n")
}

# Step 14: Check normalization information

# Search for normalization-related metadata
normalization_info <- grep(
  "normaliz|normalis|quantile|percentile|lumi",
  all_lines,
  value = TRUE,
  ignore.case = TRUE
)

# Display the beginning of matching metadata lines
for (line in normalization_info) {
  cat(substr(line, 1, 500), "\n\n")
}

# Step 15: Compare average probe expression

# Calculate mean expression for control samples
control_mean <- rowMeans(
  expression_matrix[, selected_samples$Group == "Control"]
)

# Calculate mean expression for TB samples
tb_mean <- rowMeans(
  expression_matrix[, selected_samples$Group == "TB"]
)

# Calculate the difference
mean_difference <- tb_mean - control_mean

# Summarize the differences
summary(mean_difference)

# Count positive and negative differences
table(
  Direction = ifelse(
    mean_difference > 0,
    "Higher_in_TB",
    "Lower_or_equal_in_TB"
  )
)

# Draw histogram
hist(
  mean_difference,
  breaks = 60,
  col = "lightblue",
  main = "Distribution of TB - Control Expression Differences",
  xlab = "Mean Expression Difference (TB - Control)"
)

abline(v = 0, col = "red", lwd = 2)

# Step 16: Compare sample mean expression

# Calculate mean expression for each sample
sample_means <- colMeans(expression_matrix)

# Create a data frame
mean_data <- data.frame(
  Sample_ID = colnames(expression_matrix),
  Group = selected_samples$Group,
  Mean_Expression = sample_means
)

# Summarize sample means by group
aggregate(
  Mean_Expression ~ Group,
  data = mean_data,
  FUN = summary
)

# Draw boxplot
boxplot(
  Mean_Expression ~ Group,
  data = mean_data,
  col = c("lightblue", "salmon"),
  main = "Mean Expression: TB vs Control",
  xlab = "Group",
  ylab = "Mean Expression"
)
# Step 17: Examine gender distribution

# Extract the gender metadata line
gender_line <- grep(
  '^!Sample_characteristics_ch1\\t"gender:',
  all_lines,
  value = TRUE
)

# Split into individual sample values
gender_values <- strsplit(
  gender_line,
  "\t",
  fixed = TRUE
)[[1]][-1]

# Clean the gender labels
gender_values <- gsub('^"|"$', "", gender_values)
gender_values <- sub("gender: ", "", gender_values)

# Add gender information to the full sample table
sample_info$Gender <- gender_values

# Match gender information to selected TB/Control samples
selected_samples$Gender <- sample_info$Gender[
  match(selected_samples$Sample_ID, sample_info$Sample_ID)
]

# Count gender distribution in each group
table(
  Disease = selected_samples$Group,
  Gender = selected_samples$Gender
)

# Step 18: Adjust differential expression for gender

# Convert gender into a factor
selected_samples$Gender <- factor(
  selected_samples$Gender,
  levels = c("F", "M")
)

# Create an adjusted design matrix
design_adjusted <- model.matrix(
  ~ Group + Gender,
  data = selected_samples
)

# Check design matrix
head(design_adjusted)

# Check columns
colnames(design_adjusted)

# Check whether the design matrix is full rank
qr(design_adjusted)$rank == ncol(design_adjusted)

# Fit the adjusted model
fit_adjusted <- lmFit(
  expression_matrix,
  design_adjusted
)

# Apply empirical Bayes
fit_adjusted <- eBayes(fit_adjusted)

# Summarize TB vs Control after adjusting for gender
summary(
  decideTests(
    fit_adjusted,
    coef = "GroupTB"
  )
)
# Step 19: Extract gender-adjusted DEG results

# Extract all adjusted results
deg_adjusted <- topTable(
  fit_adjusted,
  coef = "GroupTB",
  number = Inf,
  adjust.method = "BH"
)

# Display top 10 probes
head(deg_adjusted, 10)

# Count significant probes
sum(deg_adjusted$adj.P.Val < 0.05)

# Count significant probes with larger differences
sum(
  deg_adjusted$adj.P.Val < 0.05 &
    abs(deg_adjusted$logFC) > 1
)

# Summarize adjusted expression differences
summary(deg_adjusted$logFC)

# Step 20: Check ethnicity distribution

# Find ethnicity metadata
ethnicity_line <- grep(
  '^!Sample_characteristics_ch1\\t"ethnicity:',
  all_lines,
  value = TRUE
)

# Extract values
ethnicity_values <- strsplit(
  ethnicity_line,
  "\t",
  fixed = TRUE
)[[1]][-1]

# Remove quotation marks and labels
ethnicity_values <- gsub('^"|"$', "", ethnicity_values)
ethnicity_values <- sub("ethnicity: ", "", ethnicity_values)

# Add ethnicity to the full sample information
sample_info$Ethnicity <- ethnicity_values

# Match ethnicity to selected samples
selected_samples$Ethnicity <- sample_info$Ethnicity[
  match(selected_samples$Sample_ID, sample_info$Sample_ID)
]

# Display ethnicity distribution by disease group
table(
  Disease = selected_samples$Group,
  Ethnicity = selected_samples$Ethnicity
)

# Step 21: Adjust for gender and ethnicity

# Convert ethnicity to a factor
selected_samples$Ethnicity <- factor(
  selected_samples$Ethnicity
)

# Create the new design matrix
design_full <- model.matrix(
  ~ Group + Gender + Ethnicity,
  data = selected_samples
)

# Check design matrix columns
colnames(design_full)

# Check dimensions
dim(design_full)

# Check whether the design matrix is full rank
qr(design_full)$rank == ncol(design_full)

# Fit the model
fit_full <- lmFit(
  expression_matrix,
  design_full
)

# Apply empirical Bayes moderation
fit_full <- eBayes(fit_full)

# Summarize TB vs Control results
summary(
  decideTests(
    fit_full,
    coef = "GroupTB"
  )
)
# Step 22: Examine fully adjusted DEG results

# Extract TB versus Control results
deg_full <- topTable(
  fit_full,
  coef = "GroupTB",
  number = Inf,
  adjust.method = "BH"
)

# Display the top 10 probes
head(deg_full, 10)

# Summarize adjusted expression differences
summary(deg_full$logFC)

# Count probes with larger significant differences
sum(
  deg_full$adj.P.Val < 0.05 &
    abs(deg_full$logFC) > 1
)

# Count positive and negative estimated differences
table(
  Direction = ifelse(
    deg_full$logFC > 0,
    "Higher_in_TB",
    "Lower_or_equal_in_TB"
  )
)

# Step 23: Check technical metadata

# Find sample metadata field names
sample_metadata_fields <- grep(
  "^!Sample_",
  all_lines,
  value = TRUE
)

# Extract only the field names
metadata_names <- unique(
  sub("\t.*", "", sample_metadata_fields)
)

# Display available metadata fields
print(metadata_names)

# Search for possible batch-related fields
grep(
  "batch|date|scan|slide|chip|array|plate",
  metadata_names,
  value = TRUE,
  ignore.case = TRUE
)

# Step 24: Inspect technical processing information

# Extract scanning protocol
scan_lines <- grep(
  "^!Sample_scan_protocol\t",
  all_lines,
  value = TRUE
)

# Extract data processing information
processing_lines <- grep(
  "^!Sample_data_processing\t",
  all_lines,
  value = TRUE
)

# Extract hybridization protocol
hyb_lines <- grep(
  "^!Sample_hyb_protocol\t",
  all_lines,
  value = TRUE
)

# Display the first few entries
cat("\n--- SCANNING PROTOCOL ---\n")
print(substr(scan_lines[1], 1, 600))

cat("\n--- DATA PROCESSING ---\n")
print(substr(processing_lines[1], 1, 600))

cat("\n--- HYBRIDIZATION PROTOCOL ---\n")
print(substr(hyb_lines[1], 1, 600))


# Step 25: Compare technical protocols across samples

# Function to extract sample metadata values
get_metadata_values <- function(metadata_line) {
  
  values <- strsplit(
    metadata_line,
    "\t",
    fixed = TRUE
  )[[1]][-1]
  
  values <- gsub('^"|"$', "", values)
  
  return(values)
}

# Extract values for all samples
scan_values <- get_metadata_values(scan_lines[1])

processing_values <- get_metadata_values(processing_lines[1])

hyb_values <- get_metadata_values(hyb_lines[1])

# Check number of unique descriptions
cat("Unique scanning protocols:",
    length(unique(scan_values)), "\n")

cat("Unique processing descriptions:",
    length(unique(processing_values)), "\n")

cat("Unique hybridization protocols:",
    length(unique(hyb_values)), "\n")

# Check number of samples represented
cat("Number of samples:",
    length(processing_values), "\n")

# Step 26A: Load GEOquery for platform annotation

# Install GEOquery if it is not already installed
if (!requireNamespace("GEOquery", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  
  BiocManager::install("GEOquery")
}

# Load the package
library(GEOquery)

# Download platform annotation
gpl <- getGEO("GPL10558", destdir = "data")

# Extract annotation table
annotation_table <- Table(gpl)

# Check annotation dimensions
dim(annotation_table)

# Display column names
colnames(annotation_table)

# Display first few rows
head(annotation_table, 3)

# Step 26B: Match significant probes to gene symbols

# Add probe IDs to the adjusted DEG results
deg_full$Probe_ID <- rownames(deg_full)

# Select significant probes
significant_full <- deg_full[
  deg_full$adj.P.Val < 0.05 &
    abs(deg_full$logFC) > 1,
]

# Step 27: Clean gene annotations

# Keep probes with non-empty gene symbols
annotated_deg <- significant_full[
  !is.na(significant_full$Gene_Symbol) &
    trimws(significant_full$Gene_Symbol) != "",
]

# Remove extra spaces
annotated_deg$Gene_Symbol <- trimws(
  annotated_deg$Gene_Symbol
)

# Count annotated probes
nrow(annotated_deg)

# Count unique gene symbols
length(unique(annotated_deg$Gene_Symbol))

# Count genes represented by multiple probes
gene_counts <- table(annotated_deg$Gene_Symbol)

sum(gene_counts > 1)

# Display the top 10 genes with multiple probes
head(
  sort(gene_counts, decreasing = TRUE),
  10
)

# Check the number of excluded probes
nrow(significant_full) - nrow(annotated_deg)

# Check number of significant probes
nrow(significant_full)

# Match probe IDs with platform annotation
match_index <- match(
  significant_full$Probe_ID,
  annotation_table$ID
)

# Add gene symbols
significant_full$Gene_Symbol <- annotation_table$Symbol[
  match_index
]

# Add Entrez Gene IDs
significant_full$Entrez_Gene_ID <- annotation_table$Entrez_Gene_ID[
  match_index
]

# Display first 10 annotated results
head(
  significant_full[
    ,
    c(
      "Probe_ID",
      "Gene_Symbol",
      "Entrez_Gene_ID",
      "logFC",
      "adj.P.Val"
    )
  ],
  10
)

# Count probes matched to annotation records
sum(!is.na(match_index))

# Count probes with non-empty gene symbols
sum(
  !is.na(significant_full$Gene_Symbol) &
    trimws(significant_full$Gene_Symbol) != ""
)

# Step 28: Create a unique gene table

# Sort probes by adjusted p-value
annotated_deg_sorted <- annotated_deg[
  order(annotated_deg$adj.P.Val),
]

# Keep one representative probe per gene symbol
unique_genes <- annotated_deg_sorted[
  !duplicated(annotated_deg_sorted$Gene_Symbol),
]

# Check number of unique genes
nrow(unique_genes)

# Display top 10 genes
head(
  unique_genes[
    ,
    c(
      "Gene_Symbol",
      "Probe_ID",
      "Entrez_Gene_ID",
      "logFC",
      "adj.P.Val"
    )
  ],
  10
)

# Count positive and negative expression differences
table(
  Direction = ifelse(
    unique_genes$logFC > 0,
    "Higher_in_TB",
    "Lower_in_TB"
  )
)

# Save unique gene results
write.csv(
  unique_genes,
  "results/03_unique_TB_gene_results.csv",
  row.names = FALSE
)

# Step 29: Prepare background gene list

# Get all tested probe IDs
all_probe_ids <- rownames(deg_full)

# Match all tested probes to platform annotation
all_match_index <- match(
  all_probe_ids,
  annotation_table$ID
)

# Extract gene symbols
all_gene_symbols <- annotation_table$Symbol[
  all_match_index
]

# Extract Entrez Gene IDs
all_entrez_ids <- annotation_table$Entrez_Gene_ID[
  all_match_index
]

# Create background annotation table
background_annotation <- data.frame(
  Probe_ID = all_probe_ids,
  Gene_Symbol = all_gene_symbols,
  Entrez_Gene_ID = all_entrez_ids
)

# Remove missing or empty gene symbols
background_annotation <- background_annotation[
  !is.na(background_annotation$Gene_Symbol) &
    trimws(background_annotation$Gene_Symbol) != "",
]

# Count annotated background probes
nrow(background_annotation)

# Count unique background gene symbols
length(unique(background_annotation$Gene_Symbol))

# Count background entries with non-empty Entrez IDs
sum(
  !is.na(background_annotation$Entrez_Gene_ID) &
    trimws(as.character(
      background_annotation$Entrez_Gene_ID
    )) != ""
)

# Step 30: Prepare Entrez Gene IDs

# Extract Entrez IDs from significant genes
significant_entrez <- as.character(
  unique_genes$Entrez_Gene_ID
)

# Remove missing and empty IDs
significant_entrez <- significant_entrez[
  !is.na(significant_entrez) &
    trimws(significant_entrez) != ""
]

# Keep unique IDs
significant_entrez <- unique(significant_entrez)

# Extract Entrez IDs from background
background_entrez <- as.character(
  background_annotation$Entrez_Gene_ID
)

# Remove missing and empty IDs
background_entrez <- background_entrez[
  !is.na(background_entrez) &
    trimws(background_entrez) != ""
]

# Keep unique background IDs
background_entrez <- unique(background_entrez)

# Keep significant IDs present in the background
significant_entrez <- intersect(
  significant_entrez,
  background_entrez
)

# Display counts
cat(
  "Unique significant Entrez IDs:",
  length(significant_entrez),
  "\n"
)

cat(
  "Unique background Entrez IDs:",
  length(background_entrez),
  "\n"
)

# Display first 10 significant IDs
head(significant_entrez, 10)


# Step 31A: Install GO enrichment packages

# Install BiocManager if needed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Install required packages
packages <- c(
  "clusterProfiler",
  "org.Hs.eg.db",
  "enrichplot"
)

for (pkg in packages) {
  
  if (!requireNamespace(pkg, quietly = TRUE)) {
    
    BiocManager::install(
      pkg,
      ask = FALSE,
      update = FALSE
    )
    
  }
}

# Load packages
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

# Check package versions
packageVersion("clusterProfiler")
packageVersion("org.Hs.eg.db")
packageVersion("enrichplot")

# Step 31B: Validate Entrez Gene IDs

# Get all valid Entrez IDs from human annotation database
valid_entrez <- keys(
  org.Hs.eg.db,
  keytype = "ENTREZID"
)

# Keep only valid significant gene IDs
valid_significant_entrez <- intersect(
  significant_entrez,
  valid_entrez
)

# Keep only valid background gene IDs
valid_background_entrez <- intersect(
  background_entrez,
  valid_entrez
)

# Ensure significant genes belong to the background
valid_significant_entrez <- intersect(
  valid_significant_entrez,
  valid_background_entrez
)

# Display results
cat(
  "Valid significant Entrez IDs:",
  length(valid_significant_entrez),
  "\n"
)

cat(
  "Valid background Entrez IDs:",
  length(valid_background_entrez),
  "\n"
)

cat(
  "Significant IDs not recognized:",
  length(setdiff(significant_entrez, valid_entrez)),
  "\n"
)

# Step 32: GO Biological Process Enrichment Analysis

go_bp <- enrichGO(
  gene = valid_significant_entrez,
  universe = valid_background_entrez,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable = TRUE
)

# Convert enrichment results into a data frame
go_results <- as.data.frame(go_bp)

# Count significant GO biological processes
cat(
  "Number of enriched GO Biological Processes:",
  nrow(go_results),
  "\n"
)

# Display the first 10 results
head(go_results, 10)

# Step 32A: Restore significant gene results

# Load previously saved results
unique_genes <- read.csv(
  "results/03_unique_TB_gene_results.csv",
  stringsAsFactors = FALSE
)

# Check the number of genes
dim(unique_genes)

# Display the first six genes
head(unique_genes)

# Check Entrez Gene IDs
head(unique_genes$Entrez_Gene_ID)

# Step 32B: Load GEO platform annotation

library(GEOquery)

# Load GPL10558 annotation
gpl <- getGEO(
  "GPL10558",
  destdir = "data"
)

# Extract annotation table
annotation_table <- Table(gpl)

# Check dimensions
dim(annotation_table)

# Display important columns
head(
  annotation_table[
    , c("ID", "Symbol", "Entrez_Gene_ID")
  ]
)

# Step 32C: Restore expression probe IDs

# Load the GSE42830 expression matrix
expression_data <- read.delim(
  "data/GSE42830_series_matrix.txt.gz",
  comment.char = "!",
  check.names = FALSE
)

# Check dimensions
dim(expression_data)

# Display first six probe IDs
head(expression_data$ID_REF)

# Step 32D: Rebuild background Entrez Gene IDs

# Extract probe IDs from expression dataset
all_probe_ids <- expression_data$ID_REF

# Match expression probes with platform annotation
match_index <- match(
  all_probe_ids,
  annotation_table$ID
)

# Extract corresponding Entrez IDs
background_entrez <- as.character(
  annotation_table$Entrez_Gene_ID[match_index]
)

# Remove missing and empty IDs
background_entrez <- background_entrez[
  !is.na(background_entrez) &
    trimws(background_entrez) != ""
]

# Keep unique Entrez IDs
background_entrez <- unique(background_entrez)

# Extract significant Entrez IDs
significant_entrez <- as.character(
  unique_genes$Entrez_Gene_ID
)

# Remove missing and empty IDs
significant_entrez <- significant_entrez[
  !is.na(significant_entrez) &
    trimws(significant_entrez) != ""
]

# Keep unique significant IDs
significant_entrez <- unique(significant_entrez)

# Restrict significant IDs to background
significant_entrez <- intersect(
  significant_entrez,
  background_entrez
)

# Display results
cat(
  "Background Entrez IDs:",
  length(background_entrez),
  "\n"
)

cat(
  "Significant Entrez IDs:",
  length(significant_entrez),
  "\n"
)

# Step 32E: Validate Entrez Gene IDs

# Load human gene annotation database
library(org.Hs.eg.db)

# Retrieve recognized Entrez IDs
valid_entrez <- AnnotationDbi::keys(
  org.Hs.eg.db,
  keytype = "ENTREZID"
)

# Validate background Entrez IDs
valid_background_entrez <- intersect(
  background_entrez,
  valid_entrez
)

# Validate significant Entrez IDs
valid_significant_entrez <- intersect(
  significant_entrez,
  valid_background_entrez
)

# Display results
cat(
  "Valid background Entrez IDs:",
  length(valid_background_entrez),
  "\n"
)

cat(
  "Valid significant Entrez IDs:",
  length(valid_significant_entrez),
  "\n"
)

# Check unrecognized significant IDs
cat(
  "Unrecognized significant IDs:",
  length(setdiff(significant_entrez, valid_entrez)),
  "\n"
)

# Step 33: GO Biological Process Enrichment Analysis

# Load required packages
library(clusterProfiler)
library(org.Hs.eg.db)

# Perform GO Biological Process enrichment
go_bp <- clusterProfiler::enrichGO(
  gene = valid_significant_entrez,
  universe = valid_background_entrez,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable = TRUE
)

# Convert GO results into a data frame
go_results <- as.data.frame(go_bp)

# Count enriched biological processes
cat(
  "Number of enriched GO Biological Processes:",
  nrow(go_results),
  "\n"
)

# Display first 10 GO terms
head(go_results, 10)

# Step 34: GO Enrichment Dot Plot

library(enrichplot)
library(ggplot2)

# Create GO enrichment dot plot
go_dotplot <- enrichplot::dotplot(
  go_bp,
  showCategory = 15,
  x = "GeneRatio",
  color = "p.adjust",
  size = "Count",
  title = "GO Biological Process Enrichment - TB vs Control"
)

# Display plot
print(go_dotplot)

# Save plot to figures folder
dir.create("figures", showWarnings = FALSE)

ggsave(
  filename = "figures/03_GO_BP_Dotplot.png",
  plot = go_dotplot,
  width = 12,
  height = 8,
  dpi = 300
)

# Improve GO dot plot readability

library(ggplot2)

go_dotplot_final <- go_dotplot +
  scale_y_discrete(
    labels = function(x) {
      stringr::str_wrap(x, width = 35)
    }
  ) +
  theme(
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 11),
    plot.title = element_text(
      size = 14,
      face = "bold",
      hjust = 0.5
    )
  )

# Display improved plot
print(go_dotplot_final)

# Save final plot
ggsave(
  "figures/03_GO_BP_Dotplot_Final.png",
  plot = go_dotplot_final,
  width = 14,
  height = 10,
  dpi = 300
)

# Step 35: Save GO Biological Process Results

# Save all enriched GO terms
write.csv(
  go_results,
  file = "results/04_GO_BP_Enrichment_Results.csv",
  row.names = FALSE
)

# Check whether the file was saved
file.exists(
  "results/04_GO_BP_Enrichment_Results.csv"
)

# Check number of saved GO terms
nrow(go_results)


# Step 36A: KEGG Pathway Enrichment Analysis

library(clusterProfiler)

# Perform KEGG pathway enrichment
kegg_result <- enrichKEGG(
  gene = valid_significant_entrez,
  universe = valid_background_entrez,
  organism = "hsa",
  keyType = "ncbi-geneid",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.2
)

# Convert results into a data frame
kegg_results <- as.data.frame(kegg_result)

# Count enriched pathways
cat(
  "Number of enriched KEGG pathways:",
  nrow(kegg_results),
  "\n"
)

# Display first 10 pathways
head(
  kegg_results[, c(
    "ID",
    "Description",
    "GeneRatio",
    "p.adjust",
    "Count"
  )],
  10
)

# Step 36B: KEGG Pathway Dot Plot

library(enrichplot)
library(ggplot2)

# Create KEGG enrichment dot plot
kegg_dotplot <- enrichplot::dotplot(
  kegg_result,
  showCategory = 15,
  x = "GeneRatio",
  color = "p.adjust",
  size = "Count",
  title = "KEGG Pathway Enrichment - TB vs Control"
) +
  scale_y_discrete(
    labels = function(x) {
      stringr::str_wrap(x, width = 35)
    }
  ) +
  theme(
    axis.text.y = element_text(size = 10),
    plot.title = element_text(
      size = 14,
      face = "bold",
      hjust = 0.5
    )
  )

# Display plot
print(kegg_dotplot)

# Save plot
ggsave(
  "figures/04_KEGG_Pathway_Dotplot.png",
  plot = kegg_dotplot,
  width = 14,
  height = 10,
  dpi = 300
)

# Save all 44 KEGG pathway results
write.csv(
  kegg_results,
  file = "results/05_KEGG_Pathway_Results.csv",
  row.names = FALSE
)

# Verify the file
file.exists("results/05_KEGG_Pathway_Results.csv")


# Step 37A: Prepare genes for interaction network

# Load saved unique significant genes
unique_genes <- read.csv(
  "results/03_unique_TB_gene_results.csv",
  stringsAsFactors = FALSE
)

# Sort genes by adjusted p-value
sorted_genes <- unique_genes[
  order(unique_genes$adj.P.Val),
]

# Select top 30 genes
top30_genes <- head(sorted_genes, 30)

# Display selected genes
top30_genes[, c(
  "Gene_Symbol",
  "logFC",
  "adj.P.Val"
)]

# Check number of genes
nrow(top30_genes)

# Step 37B: Check STRINGdb installation

requireNamespace("STRINGdb", quietly = TRUE)

# Step 37C: Map genes to STRING protein IDs

library(STRINGdb)

# Initialize STRING database for humans
string_db <- STRINGdb$new(
  version = "12.0",
  species = 9606,
  score_threshold = 400
)

# Prepare top 30 genes
genes_for_string <- data.frame(
  Gene_Symbol = top30_genes$Gene_Symbol
)

# Map gene symbols to STRING IDs
mapped_genes <- string_db$map(
  genes_for_string,
  "Gene_Symbol",
  removeUnmappedRows = TRUE
)

# Display mapped genes
head(mapped_genes, 10)

# Count mapped genes
cat(
  "Number of genes successfully mapped:",
  nrow(mapped_genes),
  "\n"
)

# Step 37D: Generate STRING protein interaction network

# Extract mapped STRING protein IDs
protein_ids <- unique(mapped_genes$STRING_id)

# Display protein interaction network
string_db$plot_network(protein_ids)


# Step 37E: Extract STRING protein interactions

# Get interactions between mapped proteins
ppi_interactions <- string_db$get_interactions(
  protein_ids
)

# Check number of interaction records
cat(
  "Number of interaction records:",
  nrow(ppi_interactions),
  "\n"
)

# Display first 10 interactions
head(ppi_interactions, 10)

# Step 38A: Identify hub proteins

# Load igraph package
library(igraph)

# Create network from STRING interactions
ppi_graph <- graph_from_data_frame(
  d = ppi_interactions[, c("from", "to")],
  directed = FALSE,
  vertices = data.frame(
    name = protein_ids
  )
)

# Calculate degree centrality
protein_degree <- degree(ppi_graph)

# Create hub protein table
hub_results <- data.frame(
  STRING_id = names(protein_degree),
  Degree = as.numeric(protein_degree)
)

# Add gene symbols
hub_results$Gene_Symbol <- mapped_genes$Gene_Symbol[
  match(hub_results$STRING_id, mapped_genes$STRING_id)
]

# Sort by degree (highest first)
hub_results <- hub_results[
  order(-hub_results$Degree),
]

# Display top 10 hub proteins
head(
  hub_results[, c("Gene_Symbol", "Degree")],
  10
)

# Step 38B: Save hub protein results

# Save all 27 proteins and their degrees
write.csv(
  hub_results,
  "results/06_STRING_Hub_Proteins.csv",
  row.names = FALSE
)

# Select top 10 proteins
top10_hubs <- head(hub_results, 10)

# Create bar chart
library(ggplot2)

hub_plot <- ggplot(
  top10_hubs,
  aes(
    x = reorder(Gene_Symbol, Degree),
    y = Degree
  )
) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(
    title = "Top 10 Hub Proteins - TB vs Control",
    x = "Gene Symbol",
    y = "Degree Centrality"
  ) +
  theme_minimal(base_size = 13)

# Display plot
print(hub_plot)

# Save plot
ggsave(
  "figures/06_STRING_Hub_Proteins.png",
  plot = hub_plot,
  width = 10,
  height = 7,
  dpi = 300
)

# Verify saved files
file.exists("results/06_STRING_Hub_Proteins.csv")
file.exists("figures/06_STRING_Hub_Proteins.png")

# ==================================================
# Step 39A: Check Differential Expression Distribution
# ==================================================

# Load all DEG results
all_deg <- read.csv(
  "results/01_all_DEG_results.csv",
  stringsAsFactors = FALSE
)

# Count significantly upregulated probes
upregulated <- sum(
  all_deg$adj.P.Val < 0.05 &
    all_deg$logFC > 0,
  na.rm = TRUE
)

# Count significantly downregulated probes
downregulated <- sum(
  all_deg$adj.P.Val < 0.05 &
    all_deg$logFC < 0,
  na.rm = TRUE
)

# Display results
cat("Total probes:", nrow(all_deg), "\n")
cat("Upregulated probes:", upregulated, "\n")
cat("Downregulated probes:", downregulated, "\n")

# Check overall logFC distribution
summary(all_deg$logFC)


# ==================================================
# Step 39B: Check Expression Normalization
# ==================================================

# Check whether the expression matrix exists
print(exists("expression_matrix"))

# Check whether sample metadata exists
print(exists("selected_samples"))

# Check dimensions if available
if (exists("expression_matrix")) {
  print(dim(expression_matrix))
}

if (exists("selected_samples")) {
  print(dim(selected_samples))
}

# ==================================================
# Step 39C: Restore Expression Matrix
# ==================================================

# Check whether expression_data exists
print(exists("expression_data"))

# Extract expression values
expression_matrix <- as.matrix(
  expression_data[, -1]
)

# Set probe IDs as row names
rownames(expression_matrix) <- expression_data$ID_REF

# Convert values to numeric
storage.mode(expression_matrix) <- "numeric"

# Check dimensions
dim(expression_matrix)

# Display first few values
expression_matrix[1:5, 1:5]

# ==================================================
# Step 39D: Restore Sample Metadata
# ==================================================

library(GEOquery)

# Load GSE42830 from the local series matrix
gse <- getGEO(
  filename = "data/GSE42830_series_matrix.txt.gz",
  GSEMatrix = TRUE
)
# ==================================================
# Step 39E: Restore TB and Control Samples
# ==================================================

# Extract disease classifications
disease_labels <- trimws(
  as.character(sample_metadata[["disease state:ch1"]])
)

# Display all disease categories
table(disease_labels, useNA = "ifany")

# Select only TB and Control samples
selected_samples <- sample_metadata[
  !is.na(disease_labels) &
    disease_labels %in% c("TB", "Control"),
  ,
  drop = FALSE
]

# Create group labels
selected_samples$Group <- factor(
  selected_samples[["disease state:ch1"]],
  levels = c("Control", "TB")
)

# Extract corresponding expression columns
selected_expression <- expression_matrix[
  ,
  as.character(selected_samples$geo_accession),
  drop = FALSE
]

# Check sample counts
table(selected_samples$Group)

# Check expression matrix dimensions
dim(selected_expression)

# Verify sample order
identical(
  colnames(selected_expression),
  as.character(selected_samples$geo_accession)
)
# Extract sample metadata
sample_metadata <- pData(gse)

# Check dimensions
dim(sample_metadata)

# Display available metadata columns
colnames(sample_metadata)

# Display first six sample titles
head(sample_metadata$title)

# ==================================================
# Step 39F: Check Expression Distributions
# ==================================================

library(ggplot2)

# Calculate median expression for each sample
sample_medians <- apply(
  selected_expression,
  2,
  median,
  na.rm = TRUE
)

# Create data frame for visualization
median_data <- data.frame(
  Sample = names(sample_medians),
  Median_Expression = as.numeric(sample_medians),
  Group = selected_samples$Group
)

# Create plot
median_plot <- ggplot(
  median_data,
  aes(
    x = Group,
    y = Median_Expression,
    fill = Group
  )
) +
  geom_boxplot(alpha = 0.6) +
  geom_jitter(width = 0.15, size = 2) +
  labs(
    title = "Median Gene Expression: TB vs Control",
    x = "Patient Group",
    y = "Median Expression"
  ) +
  theme_minimal(base_size = 13)

# Display plot
print(median_plot)

# Display group summaries
aggregate(
  Median_Expression ~ Group,
  data = median_data,
  FUN = median
)


# ==================================================
# Step 39G: Check Original Data Preprocessing
# ==================================================

# Display preprocessing information
processing_info <- unique(
  as.character(sample_metadata$data_processing)
)

cat(
  paste(processing_info, collapse = "\n\n"),
  "\n"
)

# Display exact group median values
aggregate(
  Median_Expression ~ Group,
  data = median_data,
  FUN = median
)

# ==================================================
# Step 39H: Inspect Processing and Technical Metadata
# ==================================================

# Display processing descriptions
unique(sample_metadata$data_processing)

# Check available technical metadata columns
technical_columns <- grep(
  "batch|scan|date|platform|protocol|series",
  colnames(sample_metadata),
  value = TRUE,
  ignore.case = TRUE
)

print(technical_columns)

# Check whether processing differs between groups
table(
  selected_samples$Group,
  selected_samples$data_processing,
  useNA = "ifany"
)

# Check platform distribution
table(
  selected_samples$Group,
  selected_samples$platform_id,
  useNA = "ifany"
)

# ==================================================
# Step 39I: Check Demographic Confounding
# ==================================================

# Extract gender and ethnicity
selected_samples$Gender <- trimws(
  as.character(selected_samples[["gender:ch1"]])
)

selected_samples$Ethnicity <- trimws(
  as.character(selected_samples[["ethnicity:ch1"]])
)

# Gender distribution
cat("\nGender distribution:\n")

print(
  table(
    selected_samples$Group,
    selected_samples$Gender,
    useNA = "ifany"
  )
)

# Ethnicity distribution
cat("\nEthnicity distribution:\n")

print(
  table(
    selected_samples$Group,
    selected_samples$Ethnicity,
    useNA = "ifany"
  )
)

# Check missing values
cat(
  "\nMissing gender values:",
  sum(is.na(selected_samples$Gender) |
        selected_samples$Gender == ""),
  "\n"
)

cat(
  "Missing ethnicity values:",
  sum(is.na(selected_samples$Ethnicity) |
        selected_samples$Ethnicity == ""),
  "\n"
)

# ==================================================
# Step 39J: Validate Differential Expression Model
# ==================================================

# Set reference categories
selected_samples$Group <- factor(
  selected_samples$Group,
  levels = c("Control", "TB")
)

selected_samples$Gender <- factor(
  selected_samples$Gender
)

selected_samples$Ethnicity <- factor(
  selected_samples$Ethnicity
)

# Construct statistical design matrix
design_check <- model.matrix(
  ~ Group + Gender + Ethnicity,
  data = selected_samples
)

# Check dimensions
cat(
  "Design matrix dimensions:",
  dim(design_check),
  "\n"
)

# Check matrix rank
cat(
  "Design matrix rank:",
  qr(design_check)$rank,
  "\n"
)

# Display coefficient names
colnames(design_check)

# Check whether the model is full rank
cat(
  "Full rank:",
  qr(design_check)$rank == ncol(design_check),
  "\n"
)

# ==================================================
# Step 39K: Adjusted Sample Median Expression Analysis
# ==================================================

# Calculate sample median expression
sample_medians <- apply(
  selected_expression,
  2,
  median,
  na.rm = TRUE
)

# Create sample-level data
median_check <- data.frame(
  MedianExpression = as.numeric(sample_medians),
  Group = selected_samples$Group,
  Gender = selected_samples$Gender,
  Ethnicity = selected_samples$Ethnicity
)

# Fit adjusted linear model
median_model <- lm(
  MedianExpression ~ Group + Gender + Ethnicity,
  data = median_check
)

# Display results
summary(median_model)

# ==================================================
# Step 39L: Expression Distribution Density Plot
# ==================================================

library(ggplot2)

# Calculate expression quantiles for each sample
sample_quantiles <- apply(
  selected_expression,
  2,
  quantile,
  probs = c(0.1, 0.5, 0.9),
  na.rm = TRUE
)

# Create a summary table
quantile_data <- data.frame(
  Sample = colnames(selected_expression),
  Group = selected_samples$Group,
  Q10 = as.numeric(sample_quantiles[1, ]),
  Median = as.numeric(sample_quantiles[2, ]),
  Q90 = as.numeric(sample_quantiles[3, ])
)

# Display median quantiles by group
aggregate(
  cbind(Q10, Median, Q90) ~ Group,
  data = quantile_data,
  FUN = median
)

# ==================================================
# Step 39M: Expression Matrix Quality Check
# ==================================================

# Check missing values
cat(
  "Missing values:",
  sum(is.na(selected_expression)),
  "\n"
)

# Check infinite values
cat(
  "Infinite values:",
  sum(is.infinite(selected_expression)),
  "\n"
)

# Check overall expression range
cat(
  "Expression range:\n"
)

print(
  range(selected_expression, na.rm = TRUE)
)

# Check overall expression distribution
cat(
  "\nExpression summary:\n"
)

print(
  summary(as.vector(selected_expression))
)

# ==================================================
# Step 39N: Check Gene-Level Expression Centering
# ==================================================

# Use all 95 samples, not only TB and Control
gene_means <- rowMeans(
  expression_matrix,
  na.rm = TRUE
)

# Calculate standard deviations across samples
gene_sds <- apply(
  expression_matrix,
  1,
  sd,
  na.rm = TRUE
)

# Display gene mean distribution
cat("\nGene mean summary:\n")
print(summary(gene_means))

# Display gene standard deviation distribution
cat("\nGene standard deviation summary:\n")
print(summary(gene_sds))

# Count genes with mean approximately zero
cat(
  "\nGenes with mean close to zero:",
  sum(abs(gene_means) < 0.01),
  "\n"
)

# Count genes with standard deviation approximately one
cat(
  "Genes with SD close to one:",
  sum(abs(gene_sds - 1) < 0.05),
  "\n"
)

# ==================================================
# Step 39O: Check Sample Order and Disease Groups
# ==================================================

# Extract GSM accession numbers
sample_order_check <- data.frame(
  GSM = as.character(selected_samples$geo_accession),
  Group = selected_samples$Group
)

# Extract numeric part of GSM accession
sample_order_check$GSM_Number <- as.numeric(
  sub("GSM", "", sample_order_check$GSM)
)

# Sort by GSM accession number
sample_order_check <- sample_order_check[
  order(sample_order_check$GSM_Number),
]

# Display sample order
print(sample_order_check, row.names = FALSE)

# Count group transitions in accession order
cat(
  "\nNumber of group transitions:",
  sum(
    as.character(head(sample_order_check$Group, -1)) !=
      as.character(tail(sample_order_check$Group, -1))
  ),
  "\n"
)

# ==================================================
# Step 39P: Verify Expression and Metadata Alignment
# ==================================================

# Get sample IDs from expression matrix
expression_ids <- colnames(selected_expression)

# Get sample IDs from metadata
metadata_ids <- as.character(
  selected_samples$geo_accession
)

# Check exact sample alignment
cat(
  "Expression and metadata aligned:",
  identical(expression_ids, metadata_ids),
  "\n"
)

# Check sample counts
cat(
  "Expression samples:",
  ncol(selected_expression),
  "\n"
)

cat(
  "Metadata samples:",
  nrow(selected_samples),
  "\n"
)

# Check design matrix alignment
cat(
  "Design matrix aligned:",
  identical(
    colnames(selected_expression),
    rownames(design_check)
  ),
  "\n"
)

# Check duplicate sample IDs
cat(
  "Duplicate expression sample IDs:",
  anyDuplicated(expression_ids),
  "\n"
)

# ==================================================
# Step 39Q: Reproduce Differential Expression Analysis
# ==================================================

library(limma)

# Fit the statistical model
fit_check <- lmFit(
  selected_expression,
  design_check
)

# Apply empirical Bayes moderation
fit_check <- eBayes(fit_check)

# Extract TB vs Control results
deg_check <- topTable(
  fit_check,
  coef = "GroupTB",
  number = Inf,
  adjust.method = "BH",
  sort.by = "P"
)

# Count significant probes
cat(
  "Significant probes (FDR < 0.05):",
  sum(deg_check$adj.P.Val < 0.05),
  "\n"
)

# Count upregulated probes
cat(
  "Upregulated probes:",
  sum(
    deg_check$adj.P.Val < 0.05 &
      deg_check$logFC > 0
  ),
  "\n"
)

# Count downregulated probes
cat(
  "Downregulated probes:",
  sum(
    deg_check$adj.P.Val < 0.05 &
      deg_check$logFC < 0
  ),
  "\n"
)

# Display logFC summary
summary(deg_check$logFC)

# ==================================================
# Step 39R: Compare Saved and Recalculated DEG Results
# ==================================================

# Load previously saved results
saved_deg <- read.csv(
  "results/01_all_DEG_results.csv",
  check.names = FALSE
)

# Compare dimensions
cat("Saved DEG rows:", nrow(saved_deg), "\n")
cat("New DEG rows:", nrow(deg_check), "\n")

# Display column names
cat("\nSaved DEG columns:\n")
print(colnames(saved_deg))

cat("\nNew DEG columns:\n")
print(colnames(deg_check))

# Compare significant probe counts
cat(
  "\nSaved significant probes:",
  sum(saved_deg$adj.P.Val < 0.05, na.rm = TRUE),
  "\n"
)

cat(
  "New significant probes:",
  sum(deg_check$adj.P.Val < 0.05, na.rm = TRUE),
  "\n"
)

# ==================================================
# Step 39S: Compare Saved and Recalculated DEG Values
# ==================================================

# Add probe IDs to recalculated results
new_deg <- deg_check

new_deg$Probe_ID <- rownames(deg_check)

# Match probes using their IDs
matched_index <- match(
  saved_deg$Probe_ID,
  new_deg$Probe_ID
)

# Check matching
cat(
  "Matched probes:",
  sum(!is.na(matched_index)),
  "\n"
)

cat(
  "Unmatched probes:",
  sum(is.na(matched_index)),
  "\n"
)

# Compare logFC values
logfc_difference <- saved_deg$logFC -
  new_deg$logFC[matched_index]

cat("\nLogFC difference summary:\n")

print(summary(logfc_difference))

# Compare adjusted p-values
padj_difference <- saved_deg$adj.P.Val -
  new_deg$adj.P.Val[matched_index]

cat("\nAdjusted p-value difference summary:\n")

print(summary(padj_difference))

# Calculate correlation
cat(
  "\nLogFC correlation:",
  cor(
    saved_deg$logFC,
    new_deg$logFC[matched_index],
    use = "complete.obs"
  ),
  "\n"
)
# ==================================================
# Step 39T: Compare Differential Expression Models
# ==================================================

library(limma)

# Model 1: Group only
design_1 <- model.matrix(
  ~ Group,
  data = selected_samples
)

# Model 2: Group + Gender
design_2 <- model.matrix(
  ~ Group + Gender,
  data = selected_samples
)

# Model 3: Group + Gender + Ethnicity
design_3 <- model.matrix(
  ~ Group + Gender + Ethnicity,
  data = selected_samples
)

# Function to run limma analysis
run_deg_model <- function(design) {
  
  fit <- lmFit(selected_expression, design)
  fit <- eBayes(fit)
  
  topTable(
    fit,
    coef = "GroupTB",
    number = Inf,
    adjust.method = "BH"
  )
}

# Run three models
deg_model_1 <- run_deg_model(design_1)
deg_model_2 <- run_deg_model(design_2)
deg_model_3 <- run_deg_model(design_3)

# Compare significant probe counts
model_comparison <- data.frame(
  Model = c(
    "Group only",
    "Group + Gender",
    "Group + Gender + Ethnicity"
  ),
  
  Significant_Probes = c(
    sum(deg_model_1$adj.P.Val < 0.05),
    sum(deg_model_2$adj.P.Val < 0.05),
    sum(deg_model_3$adj.P.Val < 0.05)
  )
)

print(model_comparison)

# ==================================================
# Step 39U: Identify Original DEG Model
# ==================================================

# Function to compare saved and recalculated results
compare_deg <- function(saved, current) {
  
  # Match probe IDs
  idx <- match(
    saved$Probe_ID,
    rownames(current)
  )
  
  # Calculate differences
  logfc_diff <- saved$logFC - current$logFC[idx]
  
  # Return comparison results
  c(
    Correlation = cor(
      saved$logFC,
      current$logFC[idx]
    ),
    
    Mean_Absolute_Difference = mean(
      abs(logfc_diff)
    ),
    
    Max_Absolute_Difference = max(
      abs(logfc_diff)
    )
  )
}

# Compare all three models
comparison_results <- rbind(
  Group_Only = compare_deg(saved_deg, deg_model_1),
  
  Group_Gender = compare_deg(saved_deg, deg_model_2),
  
  Group_Gender_Ethnicity = compare_deg(saved_deg, deg_model_3)
)

print(comparison_results)


# ==================================================
# Step 39V: Adjusted Model Significant DEG Summary
# ==================================================

# Select significant probes
significant_adjusted <- deg_model_3[
  !is.na(deg_model_3$adj.P.Val) &
    deg_model_3$adj.P.Val < 0.05 &
    abs(deg_model_3$logFC) > 1,
]

# Count significant probes
cat(
  "Significant probes (FDR < 0.05 and |logFC| > 1):",
  nrow(significant_adjusted),
  "\n"
)

# Count upregulated probes
cat(
  "Upregulated:",
  sum(significant_adjusted$logFC > 1),
  "\n"
)

# Count downregulated probes
cat(
  "Downregulated:",
  sum(significant_adjusted$logFC < -1),
  "\n"
)

# Display top 10 probes
head(significant_adjusted, 10)

# ==================================================
# Step 39W: Compare Significant Probes Between Models
# ==================================================

# Original Group-only model
sig_original <- deg_model_1[
  deg_model_1$adj.P.Val < 0.05 &
    abs(deg_model_1$logFC) > 1,
]

# Adjusted model
sig_adjusted <- deg_model_3[
  deg_model_3$adj.P.Val < 0.05 &
    abs(deg_model_3$logFC) > 1,
]

# Extract probe IDs
original_ids <- rownames(sig_original)
adjusted_ids <- rownames(sig_adjusted)

# Compare probe sets
common_ids <- intersect(original_ids, adjusted_ids)

original_only <- setdiff(original_ids, adjusted_ids)
adjusted_only <- setdiff(adjusted_ids, original_ids)

# Display results
cat("Original significant probes:", length(original_ids), "\n")

cat("Adjusted significant probes:", length(adjusted_ids), "\n")

cat("Common probes:", length(common_ids), "\n")

cat("Original model only:", length(original_only), "\n")

cat("Adjusted model only:", length(adjusted_only), "\n")


# ==================================================
# Step 39X: Compare Expression Direction
# ==================================================

# Extract logFC for common probes
original_logfc <- deg_model_1[
  common_ids,
  "logFC"
]

adjusted_logfc <- deg_model_3[
  common_ids,
  "logFC"
]

# Check agreement in direction
same_direction <- sum(
  sign(original_logfc) == sign(adjusted_logfc)
)

opposite_direction <- sum(
  sign(original_logfc) != sign(adjusted_logfc)
)

# Display results
cat(
  "Common significant probes:",
  length(common_ids),
  "\n"
)

cat(
  "Same expression direction:",
  same_direction,
  "\n"
)

cat(
  "Opposite expression direction:",
  opposite_direction,
  "\n"
)

# Correlation among common probes
cat(
  "LogFC correlation:",
  cor(original_logfc, adjusted_logfc),
  "\n"
)

# ==================================================
# Step 39Y: Inspect Original Experimental Protocols
# ==================================================

# Display protocol information
protocol_columns <- c(
  "extract_protocol_ch1",
  "label_protocol_ch1",
  "hyb_protocol",
  "scan_protocol"
)

for (column in protocol_columns) {
  
  cat("\n=================================\n")
  cat("Protocol:", column, "\n")
  cat("=================================\n")
  
  print(
    unique(
      as.character(
        selected_samples[[column]]
      )
    )
  )
}

# ==================================================
# Step 39Z: Check Expression Trend by GSM Accession
# ==================================================

library(ggplot2)

# Create diagnostic data
order_qc <- data.frame(
  GSM_Number = as.numeric(
    sub("GSM", "", selected_samples$geo_accession)
  ),
  Group = selected_samples$Group,
  MedianExpression = as.numeric(
    apply(selected_expression, 2, median)
  )
)

# Plot sample median expression
ggplot(
  order_qc,
  aes(
    x = GSM_Number,
    y = MedianExpression,
    color = Group
  )
) +
  geom_point(size = 2.5) +
  geom_smooth(
    method = "lm",
    se = FALSE
  ) +
  labs(
    title = "Sample Median Expression by GSM Accession",
    x = "GSM Accession Number",
    y = "Median Gene Expression"
  ) +
  theme_minimal()

# Check correlation within each group
by(
  order_qc,
  order_qc$Group,
  function(x) {
    cor(x$GSM_Number, x$MedianExpression)
  }
)


# ==================================================
# Step 40A: Inspect GEO Series Matrix Information
# ==================================================

# Open the compressed GEO file
con <- gzfile(
  "data/GSE42830_series_matrix.txt.gz",
  open = "rt"
)

# Read the beginning of the file
file_header <- readLines(
  con,
  n = 100,
  warn = FALSE
)

close(con)

# Display relevant metadata
metadata_lines <- grep(
  "processing|normaliz|description|supplementary",
  file_header,
  value = TRUE,
  ignore.case = TRUE
)

cat(
  paste(metadata_lines, collapse = "\n"),
  "\n"
)

# ==================================================
# Step 40B: Verify GEO Expression Matrix
# ==================================================

library(GEOquery)

# Extract expression values directly from GEOquery
geo_expression <- exprs(gse)

# Check dimensions
cat(
  "GEOquery expression dimensions:",
  dim(geo_expression),
  "\n"
)

cat(
  "Manual expression dimensions:",
  dim(expression_matrix),
  "\n"
)

# Check probe IDs
cat(
  "Probe IDs identical:",
  identical(
    rownames(geo_expression),
    rownames(expression_matrix)
  ),
  "\n"
)

# Check sample IDs
cat(
  "Sample IDs identical:",
  identical(
    colnames(geo_expression),
    colnames(expression_matrix)
  ),
  "\n"
)

# Compare expression values
cat(
  "Maximum absolute difference:",
  max(
    abs(geo_expression - expression_matrix),
    na.rm = TRUE
  ),
  "\n"
)
# ==================================================
# Step 40C: Inspect GEO Metadata (Corrected)
# ==================================================

library(Biobase)

# Check object class
class(gse)

# Extract experiment information
experiment_info <- experimentData(gse)

# Display study title
cat("\nStudy Title:\n")
print(experiment_info@title)

# Display study abstract
cat("\nStudy Abstract:\n")
print(abstract(experiment_info))

# Display sample processing information
cat("\nData Processing:\n")
print(unique(sample_metadata$data_processing))


# ==================================================
# Step 40D: Save Adjusted DEG Results
# ==================================================

# Create results folder if necessary
dir.create("results", showWarnings = FALSE)

# Add probe IDs
adjusted_all <- data.frame(
  Probe_ID = rownames(deg_model_3),
  deg_model_3,
  row.names = NULL
)

# Save all adjusted DEG results
write.csv(
  adjusted_all,
  "results/07_adjusted_all_DEG_results.csv",
  row.names = FALSE
)

# Select significant adjusted probes
adjusted_significant <- subset(
  adjusted_all,
  adj.P.Val < 0.05 & abs(logFC) > 1
)

# Save significant adjusted probes
write.csv(
  adjusted_significant,
  "results/08_adjusted_significant_DEG_results.csv",
  row.names = FALSE
)

# Verify files
cat(
  "All adjusted DEG file saved:",
  file.exists("results/07_adjusted_all_DEG_results.csv"),
  "\n"
)

cat(
  "Significant adjusted DEG file saved:",
  file.exists("results/08_adjusted_significant_DEG_results.csv"),
  "\n"
)

cat(
  "Significant probes saved:",
  nrow(adjusted_significant),
  "\n"
)

# ==================================================
# Step 40E: Inspect Gene Annotation Columns
# ==================================================

# Check whether annotation table exists
cat(
  "Annotation table available:",
  exists("annotation_table"),
  "\n"
)

# Display annotation column names
if (exists("annotation_table")) {
  
  cat("\nAnnotation columns:\n")
  
  print(colnames(annotation_table))
  
  cat("\nAnnotation dimensions:\n")
  
  print(dim(annotation_table))
}


# ==================================================
# Step 40F: Annotate Adjusted Significant Probes
# ==================================================

# Select required annotation columns
gene_annotation <- annotation_table[
  ,
  c("ID", "Symbol", "Entrez_Gene_ID")
]

# Rename columns
colnames(gene_annotation) <- c(
  "Probe_ID",
  "Gene_Symbol",
  "Entrez_ID"
)

# Remove duplicate annotation probe IDs
gene_annotation <- gene_annotation[
  !duplicated(gene_annotation$Probe_ID),
]

# Match adjusted significant probes to annotation
matched_ids <- match(
  adjusted_significant$Probe_ID,
  gene_annotation$Probe_ID
)

# Combine DEG results with annotation
adjusted_annotated <- cbind(
  adjusted_significant,
  gene_annotation[
    matched_ids,
    c("Gene_Symbol", "Entrez_ID"),
    drop = FALSE
  ]
)

# Count successfully annotated probes
valid_symbols <- !is.na(adjusted_annotated$Gene_Symbol) &
  nzchar(trimws(adjusted_annotated$Gene_Symbol)) &
  adjusted_annotated$Gene_Symbol != "---"

cat(
  "Total significant probes:",
  nrow(adjusted_annotated),
  "\n"
)

cat(
  "Probes with gene symbols:",
  sum(valid_symbols),
  "\n"
)

cat(
  "Unique gene symbols:",
  length(unique(
    adjusted_annotated$Gene_Symbol[valid_symbols]
  )),
  "\n"
)

# Save annotated results
write.csv(
  adjusted_annotated,
  "results/09_adjusted_annotated_DEG_results.csv",
  row.names = FALSE
)

cat(
  "Annotated file saved:",
  file.exists(
    "results/09_adjusted_annotated_DEG_results.csv"
  ),
  "\n"
)

# ==================================================
# Step 40G: Create Unique Adjusted DEG Gene List
# ==================================================

# Keep probes with valid gene symbols
adjusted_genes <- adjusted_annotated[
  valid_symbols,
]

# Sort by adjusted p-value
adjusted_genes <- adjusted_genes[
  order(adjusted_genes$adj.P.Val),
]

# Keep one representative probe per gene
unique_adjusted_genes <- adjusted_genes[
  !duplicated(adjusted_genes$Gene_Symbol),
]

# Display results
cat(
  "Total annotated probes:",
  nrow(adjusted_genes),
  "\n"
)

cat(
  "Unique genes:",
  nrow(unique_adjusted_genes),
  "\n"
)

# Display top 10 genes
print(
  head(
    unique_adjusted_genes[
      ,
      c("Gene_Symbol", "logFC", "adj.P.Val")
    ],
    10
  )
)

# Save unique gene results
write.csv(
  unique_adjusted_genes,
  "results/10_adjusted_unique_DEG_genes.csv",
  row.names = FALSE
)

cat(
  "Unique gene file saved:",
  file.exists(
    "results/10_adjusted_unique_DEG_genes.csv"
  ),
  "\n"
)

# ==================================================
# Step 40H: Compare Original and Adjusted Gene Lists
# ==================================================

# Load original unique gene results
original_genes <- read.csv(
  "results/03_unique_TB_gene_results.csv"
)

# Display original column names
cat("\nOriginal gene file columns:\n")
print(colnames(original_genes))

# Display adjusted column names
cat("\nAdjusted gene file columns:\n")
print(colnames(unique_adjusted_genes))

# Display number of rows
cat("\nOriginal gene rows:",
    nrow(original_genes), "\n")

cat("Adjusted gene rows:",
    nrow(unique_adjusted_genes), "\n")

# ==================================================
# Step 40I: Compare Original vs Adjusted Gene Lists
# ==================================================

# Extract gene symbols
original_symbols <- unique(
  as.character(original_genes$Gene_Symbol)
)

adjusted_symbols <- unique(
  as.character(unique_adjusted_genes$Gene_Symbol)
)

# Find common genes
common_genes <- intersect(
  original_symbols,
  adjusted_symbols
)

# Find genes only in original analysis
original_only <- setdiff(
  original_symbols,
  adjusted_symbols
)

# Find genes only in adjusted analysis
adjusted_only <- setdiff(
  adjusted_symbols,
  original_symbols
)

# Display results
cat("\nOriginal unique genes:",
    length(original_symbols), "\n")

cat("Adjusted unique genes:",
    length(adjusted_symbols), "\n")

cat("Common genes:",
    length(common_genes), "\n")

cat("Only in original analysis:",
    length(original_only), "\n")

cat("Only in adjusted analysis:",
    length(adjusted_only), "\n")

# Check whether gene lists are identical
cat("Gene lists identical:",
    setequal(original_symbols, adjusted_symbols),
    "\n")


# ==================================================
# Step 40J: Compare Original and Adjusted Entrez IDs
# ==================================================

# Extract Entrez IDs
original_entrez <- unique(
  as.character(original_genes$Entrez_Gene_ID)
)

adjusted_entrez <- unique(
  as.character(unique_adjusted_genes$Entrez_ID)
)

# Remove missing and blank IDs
clean_ids <- function(x) {
  x <- trimws(x)
  x <- x[
    !is.na(x) &
      x != "" &
      x != "---"
  ]
  unique(x)
}

original_entrez <- clean_ids(original_entrez)
adjusted_entrez <- clean_ids(adjusted_entrez)

# Compare identifiers
cat("\nOriginal unique Entrez IDs:",
    length(original_entrez), "\n")

cat("Adjusted unique Entrez IDs:",
    length(adjusted_entrez), "\n")

cat("Common Entrez IDs:",
    length(intersect(
      original_entrez,
      adjusted_entrez
    )), "\n")

cat("Only in original:",
    length(setdiff(
      original_entrez,
      adjusted_entrez
    )), "\n")

cat("Only in adjusted:",
    length(setdiff(
      adjusted_entrez,
      original_entrez
    )), "\n")

cat("Entrez ID sets identical:",
    setequal(
      original_entrez,
      adjusted_entrez
    ), "\n")

# ==================================================
# Step 40K: Compare Gene Expression Changes
# ==================================================

# Match original and adjusted genes by gene symbol
comparison <- merge(
  original_genes[
    , c("Gene_Symbol", "logFC", "adj.P.Val")
  ],
  unique_adjusted_genes[
    , c("Gene_Symbol", "logFC", "adj.P.Val")
  ],
  by = "Gene_Symbol",
  suffixes = c("_Original", "_Adjusted")
)

# Calculate logFC differences
comparison$logFC_Difference <-
  comparison$logFC_Adjusted -
  comparison$logFC_Original

# Check direction changes
direction_changes <- sum(
  sign(comparison$logFC_Original) !=
    sign(comparison$logFC_Adjusted)
)

# Display results
cat("Genes compared:", nrow(comparison), "\n")

cat(
  "Correlation of logFC:",
  cor(
    comparison$logFC_Original,
    comparison$logFC_Adjusted
  ),
  "\n"
)

cat(
  "Mean absolute logFC difference:",
  mean(abs(comparison$logFC_Difference)),
  "\n"
)

cat(
  "Genes changing direction:",
  direction_changes,
  "\n"
)

# ==================================================
# Step 40L: Verify Original Gene File Model
# ==================================================

# Match original genes with adjusted genes by Probe_ID
model_check <- merge(
  original_genes[
    , c("Probe_ID", "logFC", "adj.P.Val")
  ],
  unique_adjusted_genes[
    , c("Probe_ID", "logFC", "adj.P.Val")
  ],
  by = "Probe_ID",
  suffixes = c("_Original", "_Adjusted")
)

cat(
  "Matching probe IDs:",
  nrow(model_check),
  "\n"
)

cat(
  "Maximum absolute logFC difference:",
  max(abs(
    model_check$logFC_Original -
      model_check$logFC_Adjusted
  )),
  "\n"
)

cat(
  "Maximum adjusted p-value difference:",
  max(abs(
    model_check$adj.P.Val_Original -
      model_check$adj.P.Val_Adjusted
  )),
  "\n"
)

cat(
  "All 537 probe IDs identical:",
  setequal(
    original_genes$Probe_ID,
    unique_adjusted_genes$Probe_ID
  ),
  "\n"
)

# ==================================================
# Step 40M: Verify Biological Analysis Files
# ==================================================

result_files <- c(
  "results/04_GO_BP_Enrichment_Results.csv",
  "results/05_KEGG_Pathway_Results.csv",
  "results/06_STRING_Hub_Proteins.csv"
)

# Check whether files exist
file_status <- file.exists(result_files)

# Display file status
print(
  data.frame(
    File = basename(result_files),
    Exists = file_status
  ),
  row.names = FALSE
)

# Check number of rows in available files
for (file in result_files) {
  
  if (file.exists(file)) {
    
    result <- read.csv(file)
    
    cat(
      "\n",
      basename(file),
      "- Rows:",
      nrow(result),
      "\n"
    )
  }
}

# ==================================================
# Step 40N: Inspect Biological Analysis Results
# ==================================================

# Load GO results
go_results <- read.csv(
  "results/04_GO_BP_Enrichment_Results.csv"
)

# Load KEGG results
kegg_results <- read.csv(
  "results/05_KEGG_Pathway_Results.csv"
)

# Load STRING results
string_results <- read.csv(
  "results/06_STRING_Hub_Proteins.csv"
)

# Display GO results
cat("\n========== TOP 5 GO TERMS ==========\n")

print(
  head(go_results, 5)
)

# Display KEGG results
cat("\n========== TOP 5 KEGG PATHWAYS ==========\n")

print(
  head(kegg_results, 5)
)

# Display STRING column names
cat("\n========== STRING COLUMNS ==========\n")

print(colnames(string_results))

# Display STRING results
cat("\n========== TOP 10 STRING RECORDS ==========\n")

print(
  head(string_results, 10)
)

# ==================================================
# Step 40O: Inspect Microarray Gene Background
# ==================================================

# Keep annotation for probes in the expression matrix
background_annotation <- annotation_table[
  annotation_table$ID %in% rownames(selected_expression),
]

# Extract Entrez IDs
background_entrez <- unique(
  trimws(as.character(
    background_annotation$Entrez_Gene_ID
  ))
)

# Remove missing identifiers
background_entrez <- background_entrez[
  !is.na(background_entrez) &
    background_entrez != "" &
    background_entrez != "---"
]

cat(
  "Probes in expression matrix:",
  nrow(selected_expression),
  "\n"
)

cat(
  "Probes matching GPL annotation:",
  nrow(background_annotation),
  "\n"
)

cat(
  "Unique background Entrez IDs:",
  length(background_entrez),
  "\n"
)


# ==================================================
# Step 40P: Validate Background Entrez IDs
# ==================================================

# Extract Entrez IDs from matched probes
entrez_values <- trimws(
  as.character(background_annotation$Entrez_Gene_ID)
)

# Remove missing and blank values
entrez_values <- entrez_values[
  !is.na(entrez_values) &
    entrez_values != "" &
    entrez_values != "---"
]

# Identify standard numeric Entrez IDs
numeric_entrez <- grepl("^[0-9]+$", entrez_values)

# Display results
cat(
  "Non-empty annotation entries:",
  length(entrez_values), "\n"
)

cat(
  "Standard numeric entries:",
  sum(numeric_entrez), "\n"
)

cat(
  "Non-standard entries:",
  sum(!numeric_entrez), "\n"
)

cat(
  "Unique standard Entrez IDs:",
  length(unique(entrez_values[numeric_entrez])),
  "\n"
)

# Show examples of non-standard entries
cat("\nExamples of non-standard entries:\n")

print(
  head(unique(entrez_values[!numeric_entrez]), 10)
)

# ==================================================
# Step 40Q: Create Microarray Gene Background
# ==================================================

# Extract unique numeric Entrez IDs
microarray_background <- unique(
  entrez_values[numeric_entrez]
)

# Extract significant gene Entrez IDs
significant_entrez <- unique(
  trimws(as.character(
    unique_adjusted_genes$Entrez_ID
  ))
)

# Keep valid numeric significant IDs
significant_entrez <- significant_entrez[
  !is.na(significant_entrez) &
    grepl("^[0-9]+$", significant_entrez)
]

# Check overlap
cat(
  "Microarray background genes:",
  length(microarray_background),
  "\n"
)

cat(
  "Significant Entrez genes:",
  length(significant_entrez),
  "\n"
)

cat(
  "Significant genes in background:",
  sum(significant_entrez %in% microarray_background),
  "\n"
)

cat(
  "Significant genes outside background:",
  sum(!significant_entrez %in% microarray_background),
  "\n"
)

# Save background
write.csv(
  data.frame(
    Entrez_ID = microarray_background
  ),
  "results/11_microarray_gene_background.csv",
  row.names = FALSE
)

cat(
  "Background file saved:",
  file.exists(
    "results/11_microarray_gene_background.csv"
  ),
  "\n"
)

# ==================================================
# Step 40R: GO Enrichment with Microarray Background
# ==================================================

library(clusterProfiler)
library(org.Hs.eg.db)

# Perform GO Biological Process enrichment
go_adjusted <- enrichGO(
  gene = significant_entrez,
  universe = microarray_background,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

# Convert results to a data frame
go_adjusted_results <- as.data.frame(go_adjusted)

# Display number of enriched GO terms
cat(
  "Updated enriched GO terms:",
  nrow(go_adjusted_results),
  "\n"
)

# Display top 5 GO terms
print(
  head(
    go_adjusted_results[
      ,
      c("ID", "Description", "GeneRatio",
        "BgRatio", "p.adjust")
    ],
    5
  )
)

# Save updated results
write.csv(
  go_adjusted_results,
  "results/12_GO_BP_Microarray_Background.csv",
  row.names = FALSE
)

cat(
  "Updated GO file saved:",
  file.exists(
    "results/12_GO_BP_Microarray_Background.csv"
  ),
  "\n"
)

# ==================================================
# Step 40S: Verify GO Enrichment Background
# ==================================================

# Check clusterProfiler universe setting
cat("\nEnrichment universe setting:\n")

print(
  getOption("enrichment_force_universe")
)

# Compare original and updated GO term IDs
cat(
  "\nGO term sets identical:",
  setequal(
    go_results$ID,
    go_adjusted_results$ID
  ),
  "\n"
)

# Compare adjusted p-values by GO term ID
go_comparison <- merge(
  go_results[, c("ID", "p.adjust")],
  go_adjusted_results[, c("ID", "p.adjust")],
  by = "ID",
  suffixes = c("_Original", "_Updated")
)

cat(
  "Matching GO terms:",
  nrow(go_comparison),
  "\n"
)

cat(
  "Maximum adjusted p-value difference:",
  max(abs(
    go_comparison$p.adjust_Original -
      go_comparison$p.adjust_Updated
  )),
  "\n"
)

# Inspect effective background size
cat(
  "Effective GO background:",
  unique(go_adjusted_results$UniverseSize),
  "\n"
)

# ==================================================
# Step 40T: Verify GO Background Intersection
# ==================================================

library(AnnotationDbi)
library(org.Hs.eg.db)

# Obtain genes with GO Biological Process annotations
go_bp_mapping <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = keys(org.Hs.eg.db, keytype = "ENTREZID"),
  columns = c("GO", "ONTOLOGY"),
  keytype = "ENTREZID"
)

# Keep Biological Process annotations
go_bp_mapping <- subset(
  go_bp_mapping,
  ONTOLOGY == "BP" & !is.na(GO)
)

# Unique GO-annotated Entrez IDs
go_bp_entrez <- unique(go_bp_mapping$ENTREZID)

# Intersection with microarray background
go_microarray_overlap <- intersect(
  microarray_background,
  go_bp_entrez
)

cat(
  "GO BP annotated genes in database:",
  length(go_bp_entrez),
  "\n"
)

cat(
  "Microarray background Entrez IDs:",
  length(microarray_background),
  "\n"
)

cat(
  "GO BP annotated genes on microarray:",
  length(go_microarray_overlap),
  "\n"
)

cat(
  "Effective background reported by enrichment:",
  unique(go_adjusted_results$UniverseSize),
  "\n"
)

# ==================================================
# Step 40U: KEGG Enrichment with Microarray Background
# ==================================================

library(clusterProfiler)

# Run KEGG enrichment
kegg_updated <- enrichKEGG(
  gene = significant_entrez,
  universe = microarray_background,
  organism = "hsa",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

# Convert to dataframe
kegg_updated_results <- as.data.frame(kegg_updated)

# Display results
cat(
  "Updated KEGG pathways:",
  nrow(kegg_updated_results),
  "\n"
)

# Compare original and updated pathways
cat(
  "KEGG pathway sets identical:",
  setequal(
    kegg_results$ID,
    kegg_updated_results$ID
  ),
  "\n"
)

# Check effective background
cat(
  "Updated KEGG background:",
  unique(kegg_updated_results$UniverseSize),
  "\n"
)

# Save results
write.csv(
  kegg_updated_results,
  "results/13_KEGG_Microarray_Background.csv",
  row.names = FALSE
)

cat(
  "Updated KEGG file saved:",
  file.exists(
    "results/13_KEGG_Microarray_Background.csv"
  ),
  "\n"
)

# ==================================================
# Step 40V: Save Final Research Results Summary
# ==================================================

final_summary <- data.frame(
  Measurement = c(
    "Dataset",
    "TB Samples",
    "Control Samples",
    "Total Probes",
    "Significant Probes",
    "Annotated Significant Probes",
    "Unique Gene Symbols",
    "Unique Entrez IDs",
    "GO BP Enriched Terms",
    "KEGG Enriched Pathways",
    "STRING Network Proteins"
  ),
  
  Result = c(
    "GSE42830",
    sum(selected_samples$Group == "TB"),
    sum(selected_samples$Group == "Control"),
    nrow(selected_expression),
    nrow(adjusted_significant),
    sum(valid_symbols),
    nrow(unique_adjusted_genes),
    length(significant_entrez),
    nrow(go_adjusted_results),
    nrow(kegg_updated_results),
    nrow(string_results)
  )
)

# Display summary
print(final_summary)

# Save summary
write.csv(
  final_summary,
  "results/14_Final_Research_Summary.csv",
  row.names = FALSE
)

cat(
  "Final summary saved:",
  file.exists("results/14_Final_Research_Summary.csv"),
  "\n"
)


# ==========================================
# GSE25534 Independent Validation - Step 1
# ==========================================

# Load dataset
validation_data <- read.delim(
  "data/GSE25534_intermediate_background_corrected_lowess_normalized.txt.gz",
  header = TRUE,
  check.names = FALSE
)

# Check dataset dimensions
cat("Dataset dimensions:\n")
print(dim(validation_data))

# Display first 6 column names
cat("\nColumn names:\n")
print(head(colnames(validation_data)))

# Display first 5 rows and 5 columns
cat("\nSample data:\n")
print(validation_data[1:5, 1:5])

# ==========================================
# GSE25534 Validation - Step 2
# Identify sample groups
# ==========================================

library(GEOquery)

# Load sample metadata
gse25534_geo <- getGEO(
  "GSE25534",
  GSEMatrix = TRUE,
  getGPL = FALSE,
  destdir = "data"
)

# Extract sample information
validation_metadata <- pData(gse25534_geo[[1]])

# Check metadata dimensions
print(dim(validation_metadata))

# Display sample titles and accessions
print(
  validation_metadata[
    ,
    c("geo_accession", "title")
  ],
  row.names = FALSE
)

# Find disease-related metadata columns
grep(
  "disease|status|characteristics",
  colnames(validation_metadata),
  value = TRUE,
  ignore.case = TRUE
)


# ==========================================
# GSE25534 Validation - Step 3
# Check two-channel experimental design
# ==========================================

# Check important metadata
print(
  validation_metadata[
    1:5,
    c(
      "geo_accession",
      "title",
      "status",
      "source_name_ch1",
      "source_name_ch2",
      "characteristics_ch1",
      "characteristics_ch2"
    )
  ],
  row.names = FALSE
)

# Check whether expression values represent log ratios
print(
  unique(
    validation_metadata$`data_processing`
  )
)

# ==========================================
# GSE25534 Validation - Step 4
# Prepare TB vs LTBI log ratios
# ==========================================

# Select direct TB-LTBI comparisons
direct_samples <- validation_metadata[
  (validation_metadata$source_name_ch1 == "TB" &
     validation_metadata$source_name_ch2 == "LTBI") |
    (validation_metadata$source_name_ch1 == "LTBI" &
       validation_metadata$source_name_ch2 == "TB"),
]

# Match samples to expression columns
sample_ids <- direct_samples$geo_accession

stopifnot(
  all(sample_ids %in% colnames(validation_data))
)

# Extract expression matrix
validation_expression <- as.matrix(
  validation_data[, sample_ids]
)

rownames(validation_expression) <- validation_data$ID_REF
storage.mode(validation_expression) <- "numeric"

# Show results
cat("Selected arrays:", ncol(validation_expression), "\n")

print(
  table(
    direct_samples$source_name_ch1,
    direct_samples$source_name_ch2
  )
)

cat("Expression dimensions:\n")
print(dim(validation_expression))