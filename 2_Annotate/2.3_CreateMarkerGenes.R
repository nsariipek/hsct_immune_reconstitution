# Pull the previous annotations to use as a reference in the combined analysis

# Load the libraries

library(tidyverse)
library(Seurat)

# Start with a clean slate
rm(list=ls())

# # Load your annotated Seurat object
# seu_annotated <- readRDS("/Users/dz855/Partners HealthCare Dropbox/Nurefsan Sariipek/ImmuneEscapeTP53/RDS files/seu_diet_merged.rds")

# Load the data
markers <- read.csv("/Users/dz855/Partners HealthCare Dropbox/Nurefsan Sariipek/ImmuneEscapeTP53/TP53_ImmuneEscape/2_Annotate/Markers/seu_markers_tib_top50.csv", row.names = 1)

# Remove rows where the cluster is "unidentified"
markers <- markers[markers$cluster != "Unidentified", ]
# Rename the cluster "blasts" to "undetermined"
markers$cluster <- gsub("Blasts", "Undetermined", markers$cluster)

marker_list <- split(markers$gene, markers$cluster)

# Create a data frame where each cluster has its 50 genes in a column
max_genes <- 50  # Define the number of genes per cluster
marker_df <- do.call(cbind, lapply(marker_list, function(x) {
  length(x) <- max_genes  # Pad with NA if fewer than 50 genes
  x
}))

# Assign cluster names as column headers
colnames(marker_df) <- names(marker_list)

# Save to file
write.table(marker_df, file = "/Users/dz855/Partners HealthCare Dropbox/Nurefsan Sariipek/ImmuneEscapeTP53/TP53_ImmuneEscape/2_Annotate/Markers/Cohort1_MarkerGenes.txt", sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)
