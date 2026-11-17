
# ==================================================
# TB Gene Expression Analysis
# Download Dataset: GSE42830
# ==================================================

# Create data folder
dir.create("data", showWarnings = FALSE)

# Dataset download URL
dataset_url <- paste0(
  "https://ftp.ncbi.nlm.nih.gov/geo/series/",
  "GSE42nnn/GSE42830/matrix/",
  "GSE42830_series_matrix.txt.gz"
)

# Local file path
data_file <- "data/GSE42830_series_matrix.txt.gz"

# Download dataset if it does not exist
if (!file.exists(data_file)) {
  
  download.file(
    url = dataset_url,
    destfile = data_file,
    mode = "wb",
    method = "libcurl"
  )
  
  cat("Dataset downloaded successfully!\n")
  
} else {
  
  cat("Dataset already exists.\n")
  
}

# Verify downloaded file
if (!file.exists(data_file)) {
  stop("Dataset file is missing.")
}

cat("Dataset is ready for analysis.\n")
