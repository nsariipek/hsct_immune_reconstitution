# Nurefsan Sariipek, 250127
# Merge the 80% of cells, which were excluded in 1.2_DimensionalityReduction_Clustering.R, with the annotated 20% of cells and transfer UMAP coordinates as well as celltype labels.

library(tidyverse)
library(Seurat)

# Delete environment variables
rm(list=ls())

# Load the data from the previous script
seu <- readRDS(file = "~/250108_MergedSeuratObject.rds")
seu20annotated <- readRDS(file = "~/250127_seu_annotated_merged_no-scale.rds")

# Split into 20% and 80% to make dimensionality reduction and clustering manageable
seu20annotated.cells <- colnames(seu20annotated)
seu80 <- subset(seu, cells = setdiff(colnames(seu), seu20annotated.cells))

# UMAP projection of the remaining 80% of cells (similar to https://satijalab.org/seurat/articles/integration_mapping.html#unimodal-umap-projection)
seu80 <- NormalizeData(seu80)
seu.anchors <- FindTransferAnchors(reference = seu20annotated, query = seu80, dims = 1:18, reference.reduction = "pca")
seu80 <- IntegrateEmbeddings(anchorset = seu.anchors, reference = seu20annotated, query = seu80, new.reduction.name = "ref.pca")
seu80 <- ProjectUMAP(query = seu80, query.reduction = "ref.pca", reference = seu20annotated, reference.reduction = "pca", reduction.model = "umap")

# Compare UMAPs
p1 <- DimPlot(seu20annotated, reduction = "umap", group.by = "orig.ident", shuffle = T) + theme(aspect.ratio = 1, legend.position = "none")
p2 <- DimPlot(seu80, reduction = "ref.umap", group.by = "orig.ident", shuffle = T) + theme(aspect.ratio = 1, legend.position = "none")
p1 + p2

# Transfer metadata. It's not that useful for seurat_clusters, but once we have celltype annotations, hopefully this will work well. Note: it will be useful to save the full predictions table as a csv or txt file when you do this for cell types.
predictions <- TransferData(anchorset = seu.anchors, refdata = seu_annotated$celltype, dims = 1:18)
#predictions$predicted.id <- factor(predictions$predicted.id, levels = sort(as.numeric(unique(predictions$predicted.id))))
seu80 <- AddMetaData(seu80, metadata = select(predictions, predicted.id))

# Compare UMAPs
p1 <- DimPlot(seu20annotated, reduction = "umap", group.by = "seurat_clusters", shuffle = T) + theme(aspect.ratio = 1, legend.position = "none")
p2 <- DimPlot(seu80, reduction = "ref.umap", group.by = "predicted.id", shuffle = T) + theme(aspect.ratio = 1, legend.position = "none")
p1 + p2

# Merge all data
# Check that the metadata is compatible
seu20annotated@meta.data %>% head
seu80@meta.data %>% head
seu_merge <- merge(seu20annotated, seu80)


seu_merge$UMAP_1 <- seu_merge@reductions$umap@cell.embeddings[, 1]
seu_merge$UMAP_2 <- seu_merge@reductions$umap@cell.embeddings[, 2]


