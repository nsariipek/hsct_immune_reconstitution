# Nurefsan Sariipek, updated on 250104
# Dimensionality reduction and clustering of the merged Seurat object

# Load the needed libraries
library(tidyverse)
library(Seurat)
library(readxl)
library(data.table)
library(janitor)
#library(Seurat.utils)

# Set working directory (in Terra)
setwd("/home/rstudio/TP53_ImmuneEscape/1_Seurat")

# Load the data from the previous script
seu <- readRDS(file = "~/250104_MergedSeuratObject.rds")

# Scale the data (this increases the size a lot - remove the scale.data layer before saving again)
seu <- ScaleData(seu, features = rownames(seu))
gc()
# Perform linear dimensional reduction
seu <- RunPCA(seu, features = VariableFeatures(object = seu))

# Visualize PCA results in a few different ways
DimPlot(seu, reduction = "pca") #this turned empty?
DimHeatmap(seu, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(seu, dims = 1:15, cells = 500, balanced = TRUE)

# Determine the ‘dimensionality’ of the dataset
ElbowPlot(seu)


# Find Neighbors and Cluster cells
seu <- FindNeighbors(seu, dims = 1:15)
seu <- FindClusters(seu, resolution = 1)
head(Idents(seu), 5)

# Run UMAP
seu <- RunUMAP(seu, dims = 1:30)     

# Visualize UMAP
DimPlot(seu, reduction = "umap", label = TRUE) + theme(aspect.ratio = 1)

# Visualize UMAPs with different identities
UMAP_sample <- DimPlot(seu, reduction = "umap", group.by = "orig.ident") + theme(aspect.ratio = 1)
UMAP_sample

UMAP_cohort <- DimPlot(seu, reduction = "umap", group.by = "cohort") + theme(aspect.ratio = 1)
UMAP_cohort

UMAP_status <- DimPlot(seu, reduction = "umap", group.by = "sample_status") + theme(aspect.ratio = 1)
UMAP_status

# Split the UMAPs
UMAP2 <- DimPlot(seu, reduction = "umap", group.by = "orig.ident", split.by = "cohort") + theme(aspect.ratio = 1)
UMAP2
UMAP3 <- DimPlot(seu, reduction = "umap", group.by = "orig.ident", split.by = "status") + theme(aspect.ratio = 1)
UMAP3

#To clear the memory
gc()

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


