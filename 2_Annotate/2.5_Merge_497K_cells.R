# Nurefsan Sariipek and Peter van Galen, 250127
# Merge the 80% of cells, which were excluded in 1.2_DimensionalityReduction_Clustering.R, with the annotated 20% of cells and transfer UMAP coordinates as well as celltype labels.

library(tidyverse)
library(Seurat)
library(ggsci) # for scale_color_igv
library(scattermore)

# Set working directory
setwd("~/TP53_ImmuneEscape/2_Annotate/")

# Delete environment variables
rm(list=ls())

# Load the data from the previous script
seu <- readRDS(file = "~/250108_MergedSeuratObject.rds")
seu20annotated <- readRDS(file = "~/250127_seu_annotated_merged_no-scale.rds")

# Split into 20% and 80% to make dimensionality reduction and clustering manageable
seu20annotated.cells <- colnames(seu20annotated)
seu80 <- subset(seu, cells = setdiff(colnames(seu), seu20annotated.cells))
# Check that the cells numbers make sense
identical(ncol(seu80), ncol(seu)-ncol(seu20annotated))

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
predictions <- TransferData(anchorset = seu.anchors, refdata = seu20annotated$celltype, dims = 1:18)
seu80 <- AddMetaData(seu80, metadata = select(predictions, predicted.id))

# Compare UMAPs
p1 <- DimPlot(seu20annotated, reduction = "umap", group.by = "seurat_clusters", shuffle = T) + theme(aspect.ratio = 1)
p2 <- DimPlot(seu80, reduction = "ref.umap", group.by = "predicted.id", shuffle = T) + theme(aspect.ratio = 1)
p1 + p2

# Make the metadata more compatible between objects
# Compare
seu20annotated@meta.data %>% head(n=2)
seu80@meta.data %>% head(n=2)
# Remove unnecessary data from seu20annotated
seu20annotated$RNA_snn_res.1 <- NULL
seu20annotated$seurat_clusters <- NULL
seu20annotated$cell <- NULL
# Rename seu80 celltype
colnames(seu80@meta.data) <- gsub("predicted.id", "celltype", colnames(seu80@meta.data))
# Add UMAP coordinates to metadata
seu20annotated$UMAP_1 <- seu20annotated@reductions$umap@cell.embeddings[,1]
seu20annotated$UMAP_2 <- seu20annotated@reductions$umap@cell.embeddings[,2]
seu80$UMAP_1 <- seu80@reductions$ref.umap@cell.embeddings[,1]
seu80$UMAP_2 <- seu80@reductions$ref.umap@cell.embeddings[,2]
# Now the only difference is the TNK UMAP coordinates. We don't have those for seu80, and perhaps we won't need them
setdiff(colnames(seu20annotated@meta.data), colnames(seu80@meta.data))

# Merge all cells
seu_merge <- merge(seu20annotated, seu80)
# Check that the total number of cells makes sense
ncol(seu); ncol(seu_merge)

# Set cell types as a logically ordered factor
celltypes <- c("Progenitors", "Early Erythroids", "Mid Erythroids", "Late Erythroids", "Pro Monocytes",
               "Monocytes", "Non Classical Monocytes", "cDC",  "pDC", "Pro B cells", "Pre-B", "B cells",
               "Plasma cells", "CD4 Naïve", "CD4 Effector Memory", "CD4 Memory", "Treg", "CD8 Naïve",
               "CD8 Effector", "CD8 Memory", "CD8 Exhausted",  "γδ T", "NK T", "Adaptive NK", "CD56 Bright NK",
               "CD56 Dim NK", "Cycling T-NK cells", "UD1", "UD2", "UD3")
# Check that the overlap is perfect
all(seu_merge$celltype %in% celltypes)
all(celltypes %in% seu_merge$celltype)
seu_merge$celltype <- factor(seu_merge$celltype, levels = celltypes)
seu_merge@active.ident <- seu_merge$celltype

# Join layers, which is a feature of Seurat5 I don't think we need
seu_merge <- JoinLayers(seu_merge)

# Load colors from 2.3_PvG-Colors.R
celltype_colors_df <- read.table("../celltype_colors.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE, comment.char = "")
celltype_colors <- setNames(celltype_colors_df$color, celltype_colors_df$celltype)

# Use ggplot to create a UMAP (DimPlot doesn't work b/c the reductions were lost during merge):
seu_merge@meta.data %>%
  sample_frac(1) %>%  # Randomly shuffle rows
  ggplot(aes(x = UMAP_1, y = UMAP_2, color = celltype)) +
  geom_scattermore(pointsize = 8, pixels = c(4096, 4096)) +
  scale_color_manual(values = celltype_colors) +
  theme_bw() +
  theme(aspect.ratio = 1,
        panel.grid = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 3)))


# Addition by Nurefsan, Make the T cell version 250403

Tcells <- read_rds("~/250128_Tcell_subset.rds")
Tcells@meta.data %>%
  sample_frac(1) %>%  # Randomly shuffle rows
  ggplot(aes(x = UMAP_1, y = UMAP_2, color = celltype)) +
  geom_scattermore(pointsize = 8, pixels = c(4096, 4096)) +
  scale_color_manual(values = celltype_colors) +
  theme_bw() +
  theme(aspect.ratio = 1,
        panel.grid = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 3)))


# Save as png for now. Later we should figure out how to use geom_point_rast from the library ggrastr on Terra, or save the figure on another machine
ggsave("2.5_UMAP_all_cells.pdf", width = 8, height = 4.5)

# Save to persistent disk
saveRDS(seu_merge, "~/250128_seurat_annotated_497K.rds")
#seu_merge <- readRDS("~/250127_seurat_annotated_497K.rds")

# Copy to bucket
system("gsutil cp ~/250127_seurat_annotated_497K.rds gs://fc-3783b423-62ac-4c69-8c2f-98cb0ee4503b/")



