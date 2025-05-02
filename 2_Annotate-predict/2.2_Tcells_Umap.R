# Try to re-run the UMAP for T cells
library(tidyverse)
library(dplyr)
library(Seurat)
library(BoneMarrowMap)
library(symphony)
library(RColorBrewer)
library(patchwork)
library(monocle3)

# Start with a clean slate
rm(list=ls())

# Set working directory
setwd("~/TP53_ImmuneEscape/2_Annotate-predict/")

#Load the seurat object with annotations from script 2.1
seu <- readRDS("~/250424_seu_bm_annotated.rds")

# Subset CD8 T cells
seu_subset <- subset(seu, predicted_CellType %in% c("CD8 Naive","CD8 Central Memory", "CD8 Effector Memory 1","CD8 Effector Memory 2", "CD8 Tissue Resident Memory"))

seu_subset <- NormalizeData(seu_subset)
seu_subset <- FindVariableFeatures(seu_subset)

p1 <- VariableFeaturePlot(seu_subset)
LabelPoints(plot = p1, points = head(VariableFeatures(seu_subset), 20), xnudge = 0, ynudge = 0) +
  theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5)) +
  ggtitle("Variable genes")

seu_subset <- ScaleData(seu_subset, features = rownames(seu_subset))
seu_subset <- RunPCA(seu_subset, features = VariableFeatures(object = seu_subset))


# To decide which dimension looks the best can run the followiing lines
dims_to_test <- seq(10, 20, by = 2)  # Change this as needed

plots <- list()

for (d in dims_to_test) {
  cat("Running UMAP with dims = 1:", d, "\n")
  
  # Create a copy of the object
  seu_temp <- seu_subset
  
  # Re-run neighbors and UMAP on the temp object
  seu_temp <- FindNeighbors(seu_temp, dims = 1:d, verbose = FALSE)
  seu_temp <- RunUMAP(seu_temp, reduction = "pca", dims = 1:d, return.model = TRUE, verbose = FALSE)
  
  # Create UMAP plot
  p <- DimPlot(seu_temp, reduction = "umap", group.by = "predicted_CellType", shuffle = TRUE) +
    ggtitle(paste("Dims: 1-", d)) +
    theme(aspect.ratio = 1)
  
  plots[[paste0("dims_", d)]] <- p
}


if (!dir.exists("umap_dim_checks")) dir.create("umap_dim_checks")

# Save each plot
for (name in names(plots)) {
  ggsave(
    filename = paste0("umap_dim_checks/", name, ".pdf"),
    plot = plots[[name]],
    width = 6,
    height = 6
  )
}

# When you decide on the dimension run these
# Cluster data
seu_subset <- FindNeighbors(seu_subset, dims = 1:18)

seu_subset <- RunUMAP(seu_subset, reduction = "pca", dims = 1:18, return.model = T)

# Skip
#seu_subset <- FindClusters(seu_subset, resolution = 1)

# UMAP
DimPlot(seu_subset, reduction = "umap", group.by = "predicted_CellType", shuffle = T) +
  theme(aspect.ratio = 1)


