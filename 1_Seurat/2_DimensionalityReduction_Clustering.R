# Nurefsan Sariipek, updated on 250104
# Dimensionality reduction and clustering of the merged Seurat object

# Load the needed libraries
library(tidyverse)
library(Seurat)
library(readxl)
library(data.table)
library(janitor)
library(future.apply)

# Start with a clean slate
rm(list=ls())

# Set working directory (in Terra)
setwd("/home/rstudio/TP53_ImmuneEscape/1_Seurat")

# Load the data from the previous script
seu <- readRDS(file = "~/250108_MergedSeuratObject.rds")

#Since the regular downstream analysis failed many times Nurefsan tried this tutorial from Seurat:https://satijalab.org/seurat/archive/v4.3/integration_large_datasets
# First perform standard normalization and variable feature selection
seu.list <- SplitObject(seu, split.by = "orig.ident")

# Use multiple cores
plan(multisession)
# This only needs to be run once (but it doesn't hurt to do it again)
options(future.globals.maxSize = 160 * 1024^3) # 160 GB

seu.list <- future_lapply(X = seu.list, FUN = function(x) {
  print(as.character(x@meta.data$orig.ident[1]))
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, verbose = FALSE)
})

# Select features for downstream integration, and run PCA on each object in the list, which is required for running the alternative reciprocal PCA workflow

features <- SelectIntegrationFeatures(object.list = seu.list)
seu.list <- future_lapply(X = seu.list, FUN = function(x) {
  print(as.character(x@meta.data$orig.ident[1]))
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
},future.seed = T)

# Revert processing to sequential execution
plan(sequential)

anchors <- FindIntegrationAnchors(object.list = seu.list, reduction = "rpca", dims = 1:50)
seu.integrated <- IntegrateData(anchorset = anchors, dims = 1:50)
seu.integrated <- ScaleData(seu.integrated, verbose = FALSE)
seu.integrated <- RunPCA(seu.integrated, verbose = FALSE)
seu.integrated <- RunUMAP(seu.integrated, dims = 1:50)

# This created an object called seu.integrated(46.4 GB) and Nurefsan saved this 
saveRDS(seu.integrated, file = "~/250106_IntegratedSeuratObject.rds")

# We should save the Seurat object at this point (saveRDS), but first reduce the size
# @Nurefsan: I have not tested these following methods (you should)
seu_diet <- DietSeurat(seu, dimreducs = names(seu@reductions))
# Also remove the scale.data layer, which takes a lot of space
seu_filtered <- removeLayersByPattern(seu_diet, pattern = "scale.data", perl = TRUE)
# Save for quick loading later
saveRDS(seu_filtered, file = "~/250105_ClusteredSeuratObject.rds")
# @Nurefsan in downstream analyses, you should also save 2501xx_AnnotatedSeuratObject.rds, 2501xx_SeuratObject_SouporCell.rds, 2501xx_SeuratObject_numbat.rds, etc.

# Run Clusters to do cell annotation(this is the most important step and it takes some time+memory)
seu_markers <- FindAllMarkers(seu, min.pct = .3, logfc.threshold = .3)

# Save as tibble for next time
# @Nurefsan this code saves the markers as a csv file, not a tibble (which is fine, but you should update the explanation).
seu_markers_tib <- as_tibble(seu_markers)
write.csv(seu_markers_tib, file = "seu_markers_tib.csv")
# @Nurefsan I think I would prefer the following, but would have to check output to make sure
write.tsv(as_tibble(seu_markers), file = "marker_genes.tsv")


