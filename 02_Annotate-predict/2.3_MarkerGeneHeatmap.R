# Peter van Galen, 250902
# Plot a heatmap of cell type-specific marker genes

# Load libraries
library(Seurat)
library(tidyverse)

# Set working directory
# fmt: skip
setwd("~/DropboxMGB/Projects/ImmuneEscapeTP53/hsct_immune_reconstitution/10_Miscellaneous")

# Delete environment variables
rm(list = ls())

# Load data
seu <- readRDS("../AuxiliaryFiles/250528_Seurat_complete.rds")

# Subset for cells with successful annotation
seu_subset <- subset(seu, !is.na(celltype))

# Find markers for all cell types
Idents(seu_subset) <- "celltype"
seu_subset <- NormalizeData(seu_subset)
markers <- FindAllMarkers(
  seu_subset,
  only.pos = TRUE,
  min.pct = 0.1,
  logfc.threshold = 0.1
)


markers_tib <- as_tibble(markerGenes, rownames = "gene") %>%
  arrange(-avg_log2FC)
write_tsv(markers_tib, file = "2.3_Celltype_markers.txt")


markers_tib %>% pivot_wider(id_cols = log2FC, names_from, )

# Define signatures
bpdcn_sign <- markers_tib %>%
  filter(p_val_adj < 1E-30, avg_log2FC > 1) %>%
  .$gene


# Get top markers per cell type
top_markers <- markers %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 5) %>%
  pull(gene) %>%
  unique()

# Create heatmap
DoHeatmap(seu, features = top_markers, group.by = "celltype") +
  theme(axis.text.y = element_text(size = 8))
