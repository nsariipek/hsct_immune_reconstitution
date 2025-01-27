# Nurefsan Sariipek, updated on 250104
# Dimensionality reduction and clustering of the merged Seurat object

# Load the needed libraries
library(tidyverse)
library(Seurat)
#library(readxl)
#library(data.table)
#library(janitor)

# Start with a clean slate
rm(list=ls())

# Set working directory (in Terra)
setwd("/home/rstudio/TP53_ImmuneEscape/1_Seurat")
# Google VM:
#setwd("/home/unix/vangalen/TP53_ImmuneEscape/1_Seurat")

# Load the data from the previous script
seu <- readRDS(file = "~/250108_MergedSeuratObject.rds")

# Split into 20% and 80% to make dimensionality reduction and clustering manageable
seu20.cells <- colnames(seu)[seq(1, length(colnames(seu)), by = 5)]
seu80.cells <- setdiff(colnames(seu), seu20.cells)
seu20 <- subset(seu, cells = seu20.cells)
seu80 <- subset(seu, cells = seu80.cells)

# Find variable features on the smaller object
seu20 <- NormalizeData(seu20)
seu20 <- FindVariableFeatures(seu20)
p1 <- VariableFeaturePlot(seu20)
LabelPoints(plot = p1, points = head(VariableFeatures(seu20), 20), xnudge = 0, ynudge = 0) +
  theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5)) +
  ggtitle("Variable genes")

# Regular dimensionality reduction. Nurefsan: you can try different `dims` in RunUMAP to see what you like.
seu20 <- ScaleData(seu20, features = rownames(seu20))
seu20 <- RunPCA(seu20, features = VariableFeatures(object = seu20))
ElbowPlot(seu20)
# Assign number to n_dim because it makes sense to use the same number of dimensions for various downstream steps: dimensionality reduction, clustering, and integration. return.model = T is required for integration later
n_dim <- 18
seu20 <- RunUMAP(seu20, reduction = "pca", dims = 1:n_dim, return.model = T)

# Cluster data
seu20 <- FindNeighbors(seu20, dims = 1:n_dim)
seu20 <- FindClusters(seu20, resolution = 1)

# Visualize UMAPs with different groupings
DimPlot(seu20, reduction = "umap", group.by = "orig.ident", shuffle = T) + theme(aspect.ratio = 1, legend.position = "none")
DimPlot(seu20, reduction = "umap", group.by = "cohort", shuffle = T) + theme(aspect.ratio = 1)
DimPlot(seu20, reduction = "umap", group.by = "sample_status", shuffle = T) + theme(aspect.ratio = 1)
DimPlot(seu20, reduction = "umap", group.by = "seurat_clusters", shuffle = T) + theme(aspect.ratio = 1)
FeaturePlot(seu20, features = "CD34") + theme(aspect.ratio = 1)

# Save seu20 to work on annotations
saveRDS(seu20, "~/250113_SplittedSeuratObject.rds")



