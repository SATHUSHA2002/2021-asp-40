
# ============================================================
# TB Gene Expression Analysis
# Script 08: STRING Protein Interaction Network
# ============================================================

source("code/05_gene_annotation.R")

library(STRINGdb)
library(igraph)
library(ggplot2)

# 1. Select top 30 significant genes
top_genes <- head(
  unique_deg[order(unique_deg$adj.P.Val), ],
  30
)

# 2. Initialize STRING database
string_db <- STRINGdb$new(
  version = "12.0",
  species = 9606,
  score_threshold = 400
)

# 3. Map gene symbols to STRING protein IDs
mapped_genes <- string_db$map(
  top_genes,
  "Gene_Symbol",
  removeUnmappedRows = TRUE
)

cat("\nMapped proteins:",
    nrow(mapped_genes), "\n")

# 4. Get interactions
interactions <- string_db$get_interactions(
  mapped_genes$STRING_id
)

# Keep interactions among selected proteins
selected_ids <- unique(mapped_genes$STRING_id)

interactions <- interactions[
  interactions$from %in% selected_ids &
    interactions$to %in% selected_ids,
  ,
  drop = FALSE
]

# 5. Create network
network <- graph_from_data_frame(
  interactions[, c("from", "to")],
  directed = FALSE,
  vertices = data.frame(name = selected_ids)
)

network <- simplify(
  network,
  remove.multiple = TRUE,
  remove.loops = TRUE
)

# 6. Calculate hub proteins
protein_degree <- degree(network)

hub_table <- data.frame(
  STRING_id = names(protein_degree),
  Degree = as.integer(protein_degree)
)

hub_table <- merge(
  hub_table,
  mapped_genes[, c("STRING_id", "Gene_Symbol")],
  by = "STRING_id",
  all.x = TRUE
)

hub_table <- hub_table[
  order(-hub_table$Degree),
  ,
  drop = FALSE
]

# 7. Save results
write.csv(
  hub_table,
  "results/08_STRING_hub_proteins.csv",
  row.names = FALSE
)

# 8. Display summary
cat("\n===== STRING NETWORK SUMMARY =====\n")

cat("Proteins:", vcount(network), "\n")
cat("Interactions:", ecount(network), "\n")

cat("\nTop 10 hub proteins:\n")

print(
  head(
    hub_table[, c("Gene_Symbol", "Degree")],
    10
  )
)

# 9. Save network visualization
png(
  "figures/08_STRING_network.png",
  width = 1600,
  height = 1200,
  res = 150
)

set.seed(123)

plot(
  network,
  layout = layout_with_fr(network),
  vertex.label = hub_table$Gene_Symbol[
    match(V(network)$name, hub_table$STRING_id)
  ],
  vertex.size = 8 + degree(network),
  vertex.label.cex = 0.7,
  vertex.color = "skyblue",
  edge.color = "gray70",
  main = "STRING Protein Interaction Network - TB"
)

dev.off()

cat("\nSTRING network analysis completed!\n")
