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

# UMAP projection of the remaining 80% of cells (similar to https://satijalab.org/seurat/articles/integration_mapping.html#unimodal-umap-projection)
seu80 <- NormalizeData(seu80)
seu.anchors <- FindTransferAnchors(reference = seu20, query = seu80, dims = 1:n_dim, reference.reduction = "pca")
seu80 <- IntegrateEmbeddings(anchorset = seu.anchors, reference = seu20, query = seu80, new.reduction.name = "ref.pca")
seu80 <- ProjectUMAP(query = seu80, query.reduction = "ref.pca", reference = seu20, reference.reduction = "pca", reduction.model = "umap")

# Compare UMAPs
p1 <- DimPlot(seu20, reduction = "umap", group.by = "orig.ident", shuffle = T) + theme(aspect.ratio = 1, legend.position = "none")
p2 <- DimPlot(seu80, reduction = "ref.umap", group.by = "orig.ident", shuffle = T) + theme(aspect.ratio = 1, legend.position = "none")
p1 + p2

# Transfer metadata. It's not that useful for seurat_clusters, but once we have celltype annotations, hopefully this will work well. Note: it will be useful to save the full predictions table as a csv or txt file when you do this for cell types.
predictions <- TransferData(anchorset = seu.anchors, refdata = seu20$seurat_clusters, dims = 1:n_dim)
predictions$predicted.id <- factor(predictions$predicted.id, levels = sort(as.numeric(unique(predictions$predicted.id))))
seu80 <- AddMetaData(seu80, metadata = select(predictions, predicted.id))

# Compare UMAPs
p1 <- DimPlot(seu20, reduction = "umap", group.by = "seurat_clusters", shuffle = T) + theme(aspect.ratio = 1, legend.position = "none")
p2 <- DimPlot(seu80, reduction = "ref.umap", group.by = "predicted.id", shuffle = T) + theme(aspect.ratio = 1, legend.position = "none")
p1 + p2

# Once the UMAP coordinates and cell type annotations are finalized, we can merge the seu20 and seu80 (I recommend storing the UMAP coordinates in new metadata colums UMAP_1 and UMAP2 and the cell type annotations in another metadata column...can help with this)














# THE FOLLOWING IS OBSOLETE

# Nurefsan and Peter tried this tutorial from Seurat (https://satijalab.org/seurat/archive/v4.3/integration_large_datasets), but this ultimately failed
# First perform standard normalization and variable feature selection
seu.list <- SplitObject(seu, split.by = "orig.ident")

# Use multiple cores
library(future.apply)
plan(multisession)
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
gc()

# The following took Peter ~7 hours (80 vCPU, 40 core, 320 GB memory). Memory and storage can be monitored in terminal using `free -h` or `htop`. The only time we got something to work, is when we added the argument `reference = c(1,2)` but there is no good justification for that approach
anchors <- FindIntegrationAnchors(object.list = seu.list, reduction = "rpca", dims = 1:50)
# Save result in the VM and Terra's bucket (delete later to save space)
saveRDS(anchors, file = "~/250111_Anchors.rds")
system("sudo gsutil cp ~/250111_Anchors.rds gs://fc-3783b423-62ac-4c69-8c2f-98cb0ee4503b/")

# Continue with IntegrateData This yields the warnings "Warning: Layer counts isn't present in the assay object; returning NULL" and "Warning: Different cells in new layer data than already exists for scale.data". Eventually, this crashed even with 520 GB memory
seu.integrated <- IntegrateData(anchorset = anchors, dims = 1:50)
seu.integrated <- ScaleData(seu.integrated, verbose = FALSE)
seu.integrated <- RunPCA(seu.integrated, verbose = FALSE)
seu.integrated <- RunUMAP(seu.integrated, dims = 1:50)

# The seu.integrated saved here was generated with `reference = c(1,2)` (FindIntegrationAnchors) so we should not move forward with it
saveRDS(seu.integrated, file = "~/250106_IntegratedSeuratObject.rds")

# Peter is unsure if we ran the following
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


