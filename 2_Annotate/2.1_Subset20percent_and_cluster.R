# Nurefsan Sariipek, updated on 250104
# Dimensionality reduction and clustering of the merged Seurat object

# Load the needed libraries
library(tidyverse)
library(Seurat)

# Start with a clean slate
rm(list=ls())

# Set working directory (in Terra)
setwd("/home/rstudio/TP53_ImmuneEscape/2_Annotate")

# Load the data from 1.1_CreateSeuratObject.R
seu <- readRDS(file = "~/250409_MergedSeuratObject.rds")

# Select 20% of the cells to make dimensionality reduction and clustering manageable. Later on, the remaining 80% of cells are added back.
# An alternative method is described here: https://satijalab.org/seurat/articles/seurat5_sketch_analysis.html
seu20.cells <- colnames(seu)[seq(1, length(colnames(seu)), by = 5)]
seu20 <- subset(seu, cells = seu20.cells)

# Find variable features on the smaller object
seu20 <- NormalizeData(seu20)
seu20 <- FindVariableFeatures(seu20)
p1 <- VariableFeaturePlot(seu20)
LabelPoints(plot = p1, points = head(VariableFeatures(seu20), 20), xnudge = 0, ynudge = 0) +
  theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5)) +
  ggtitle("Variable genes")

# Regular dimensionality reduction. At least 130 GB memory is required.
seu20 <- ScaleData(seu20, features = rownames(seu20))
seu20 <- RunPCA(seu20, features = VariableFeatures(object = seu20))
ElbowPlot(seu20)

# Calculate UMAP coordinates. return.model = T is required for integration later
seu20 <- RunUMAP(seu20, reduction = "pca", dims = 1:18, return.model = T)

# Cluster data
seu20 <- FindNeighbors(seu20, dims = 1:18)
seu20 <- FindClusters(seu20, resolution = 1)

# Visualize UMAPs with different groupings
DimPlot(seu20, reduction = "umap", group.by = "patient_id", shuffle = T) +
  theme(aspect.ratio = 1, legend.position = "none")
DimPlot(seu20, reduction = "umap", group.by = "cohort", shuffle = T) + theme(aspect.ratio = 1)
DimPlot(seu20, reduction = "umap", group.by = "sample_status", shuffle = T) + theme(aspect.ratio = 1)
DimPlot(seu20, reduction = "umap", group.by = "seurat_clusters", shuffle = T) + theme(aspect.ratio = 1)
FeaturePlot(seu20, features = "CD34") + theme(aspect.ratio = 1)

# Save seu20 to work on annotations. Save storage by removing the scale data slot
seu20@assays$RNA$scale.data <- NULL
saveRDS(seu20, "~/250410_SubsettedSeuratObject.rds")


