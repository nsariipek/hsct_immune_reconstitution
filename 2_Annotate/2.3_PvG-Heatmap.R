# Peter van Galen, 250125
# Visualize marker expression of T cell subsets in a heatmap
# If helpful, this code should be integrated with 2.2_Tcell_Annotation.R

library(tidyverse)
library(Seurat) # Using 5.1.0
library(data.table)
library(ggpubr)

# Set working directory
setwd("~/TP53_ImmuneEscape/2_Annotate/")

# Clear environment variables
rm(list=ls())

# Store merged object locally (this is a temporary solution while analysis is ongoing - eventually we should probably put finished objects on Figshare or publish with Terra). Terminal:
#cd ~
#gsutil cp gs://fc-3783b423-62ac-4c69-8c2f-98cb0ee4503b/250124_seu_annotated_merged.rds .
# When completing analysis, delete file from local home directory (e.g. in R, `unlink("~/250124_seu_annotated_merged.rds")`)

# Load Seurat object into R memory
seu <- readRDS("~/250124_seu_annotated_merged.rds")

# Make the cell type metadata a logically ordered factor
seu$celltype %>% unique
celltypes <- c("Progenitors", "Early Erythroids", "Mid Erythroids", "Late Erythroids", "Pro Monocytes", "Monocytes", "Non Classical Monocytes", "cDC",  "pDC", "Pro B cells", "Pre-B", "B cells", "Plasma cells", "CD4 Naïve", "CD4 Effector Memory", "CD4 Memory", "Treg", "CD 8 Naïve", "CD8 Effector", "CD8 Memory", "CD8 Exhausted",  "γδ T", "NK T", "Adaptive NK", "CD56 Bright NK", "CD56 Dim NK", "Cycling T-NK cells", "UD1", "UD2", "UD3")
# Check & replace celltype metadata
all(seu$celltype %in% celltypes)
all(celltypes %in% seu$celltype)
seu$celltype <- factor(seu$celltype, levels = celltypes)

# Subset for TNK cells
TNK_celltypes <- c("CD4 Naïve", "CD4 Effector Memory", "CD4 Memory", "Treg", "CD 8 Naïve", "CD8 Effector", "CD8 Memory", "CD8 Exhausted",  "γδ T", "NK T", "Adaptive NK", "CD56 Bright NK", "CD56 Dim NK")
seuTNK <- subset(seu, celltype %in% TNK_celltypes)
# Check
DimPlot(seuTNK)

# Define marker genes
seu@active.ident <- seu$celltype
markerGenes <- FindAllMarkers(seuTNK, only.pos = T)
# Select top n genes for each cluster
top_n_tib <- as_tibble(markerGenes) %>% group_by(cluster) %>%
  arrange(-avg_log2FC, .by_group = T) %>% slice_max(order_by = avg_log2FC, n = 8)

# Plot heatmap
DoHeatmap(seuTNK, features = unique(top_n_tib$gene), draw.lines = FALSE) +
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 5, name = "RdBu")))

#ggsave("tmp.png", width = 20, height = 12)

# Delete the annotated merged Seurat object
#unlink("~/250124_seu_annotated_merged.rds")