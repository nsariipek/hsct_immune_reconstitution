# Nurefsan Sariipek and Peter van Galen, 250127
# Merge the 80% of cells, which were excluded in 1.2_DimensionalityReduction_Clustering.R, with the annotated 20% of cells and transfer UMAP coordinates as well as celltype labels.

library(tidyverse)
library(Seurat)
#library(ggsci) # for scale_color_igv
#library(scattermore)

# Set working directory
setwd("~/TP53_ImmuneEscape/2_Annotate/")

# Delete environment variables
rm(list=ls())

# Load the data from previous scripts
seu_merge <- readRDS(file = "~/250409_MergedSeuratObject.rds")
seu20_anno <- readRDS(file = "~/250411_SubsettedAnnotatedSeuratObject.rds")

# Remove the annotated cells from the full object
seu20_anno.cells <- colnames(seu20_anno)
seu80 <- subset(seu_merge, cells = setdiff(colnames(seu_merge), seu20_anno.cells))
# Check that the cells numbers make sense
identical(ncol(seu80), ncol(seu_merge)-ncol(seu20_anno))

# PROJECT FULL UMAP AND PREDICT CELL TYPES -------------------------------------
# UMAP projection of the remaining 80% of cells (similar to https://satijalab.org/seurat/articles/integration_mapping.html#unimodal-umap-projection)
seu80 <- NormalizeData(seu80)
seu.anchors <- FindTransferAnchors(reference = seu20_anno, query = seu80, dims = 1:18, reference.reduction = "pca")
seu80 <- IntegrateEmbeddings(anchorset = seu.anchors, reference = seu20_anno, query = seu80, new.reduction.name = "ref.pca")
seu80 <- ProjectUMAP(query = seu80, query.reduction = "ref.pca", reference = seu20_anno, reference.reduction = "pca", reduction.model = "umap", reduction.name = "ref.umap", reduction.key = "refUMAP_")

# Transfer celltype annotations
predictions_df <- TransferData(anchorset = seu.anchors, refdata = seu20_anno$celltype, dims = 1:18)
colnames(predictions_df) <- gsub("predicted.id", "celltype", colnames(predictions_df))
predictions_df$cell <- rownames(predictions_df)
write_csv(predictions_df[,c("cell", "celltype")], file = "2.3_CellTypePredictions.csv")
seu80 <- AddMetaData(seu80, metadata = select(predictions_df, celltype))

# Compare UMAPs
p1 <- DimPlot(seu20_anno, reduction = "umap", group.by = "celltype", shuffle = T) + theme(aspect.ratio = 1)
p2 <- DimPlot(seu80, reduction = "ref.umap", group.by = "celltype", shuffle = T) + theme(aspect.ratio = 1)
p1 + p2

# Merge all cells without dimensionality reductions (but with cell type annotations)
seu_merge_anno <- merge(seu20_anno, seu80, merge.dr = F)

# Set cell types as a logically ordered factor
celltypes <- c("Progenitors", "Early Erythroids", "Mid Erythroids", "Late Erythroids", "Pro Monocytes",
               "Monocytes", "Non Classical Monocytes", "cDC",  "pDC", "Pro B cells", "Pre-B", "B cells",
               "Plasma cells", "CD4 Naïve", "CD4 Memory", "CD4 Effector Memory", "Treg", "CD8 Naïve",
               "CD8 Effector", "CD8 Memory", "CD8 Exhausted",  "γδ T", "NK T", "Adaptive NK", "CD56 Bright NK",
               "CD56 Dim NK", "Cycling T-NK cells", "UD1", "UD2", "UD3")
# Check that the overlap is perfect
all(seu_merge_anno$celltype %in% celltypes)
all(celltypes %in% seu_merge_anno$celltype)
seu_merge_anno$celltype <- factor(seu_merge_anno$celltype, levels = celltypes)
seu_merge_anno@active.ident <- seu_merge_anno$celltype

# Layers are a Seurat v5 feature we don't need
seu_merge_anno <- JoinLayers(seu_merge_anno)

# Check that the total number of cells makes sense
ncol(seu_merge); ncol(seu_merge_anno)

# Used below
tnk_celltypes <- c("CD4 Memory", "CD8 Memory", "CD4 Naïve", "Treg", "CD8 Effector", "CD4 Effector Memory", "γδ T", "CD8 Exhausted", "CD56 Dim NK", "CD8 Naïve", "NK T", "CD56 Dim NK", "CD56 Bright NK", "Adaptive NK")
tnk_cells <- colnames( subset(seu_merge_anno, celltype %in% tnk_celltypes) )


# Project TNK UMAP -------------------------------------------------------------
# First I need to reload & subset the original merged Seurat object
seu_for_TNK <- readRDS(file = "~/250409_MergedSeuratObject.rds")
seu80_TNK <- subset(seu_for_TNK, cells = intersect(setdiff(colnames(seu_for_TNK), colnames(seu20_anno.cells)), tnk_cells))
seu80_TNK <- NormalizeData(seu80_TNK)
# I also need the TNK object saved in 2.2_Tcell_Annotation.R
seu20_T <- readRDS("250411_AnnotatedSeurat_TNK.rds")
# Then predict UMAP
seu20_T <- ScaleData(seu20_T) # not sure why, but the next line gives an error without this

seu20_T@reductions$umap <- NULL
seu20_T <- RunUMAP(seu20_T, dims = 1:20, return.model = T)
seu20_T@assays$RNA@layers$scale.data <- NULL # No idea why this would work --- I'm desparate

seuT.anchors <- FindTransferAnchors(reference = seu20_T, query = seu80_TNK, dims = 1:20, reference.reduction = "pca")
seu80_TNK <- IntegrateEmbeddings(anchorset = seuT.anchors, reference = seu20_T, query = seu80_TNK, new.reduction.name = "ref.pca")
seu80_TNK <- ProjectUMAP(query = seu80_TNK, query.reduction = "ref.pca", reference = seu20_T, reference.reduction = "pca", reduction.model = "umap", reduction.name = "ref.umapTNK", reduction.key = "refUMAPTNK_")
DimPlot(seu80_TNK, reduction = "ref.umapTNK", group.by = "patient_id") + theme(aspect.ratio = 1)


# The following also does not work because seu20_anno does not have a model saved with it's TNK UMAP
seu80_TNK2 <- subset(seu80, celltype %in% tnk_celltypes)
seu80_TNK2 <- ProjectUMAP(query = seu80_TNK2, query.reduction = "ref.pca", reference = seu20_anno, reference.reduction = "pca", reduction.model = "umapTNK", reduction.name = "ref.umapTNK", reduction.key = "refUMAPTNK_")


#seu80_tnk <- subset(seu80, celltype %in% tnk_celltypes)
#seu <- readRDS("~/250410_SubsettedSeuratObject.rds")
#seu20_anno_tnk <- subset(seu20_anno, celltype %in% tnk_celltypes)

#seu20_anno_tnk2 <- seu20_anno_tnk
#seu20_anno_tnk2 <- subset(seu20_anno, cells = seu20_tnk_cells)
#seu20_anno_tnk2 <- FindVariableFeatures(seu20_anno_tnk2)
#seu20_anno_tnk2@reductions$umapTNK <- NULL
#seu20_anno_tnk2 <- FindNeighbors(seu20_anno_tnk2, dims = 1:20)
#seu20_anno_tnk2 <- RunUMAP(seu20_anno_tnk2, dims = 1:20, reduction.name = "umapTNK")

#seu20_anno_tnk@reductions
#DimPlot(seu20_anno_tnk, reduction = "umapTNK")
#DimPlot(seu20_anno_tnk2, reduction = "umapTNK")

#seu.anchors <- FindTransferAnchors(reference = seu20_anno_tnk, query = seu80_tnk, dims = 1:20, reference.reduction = "pca")
#seu80_tnk <- IntegrateEmbeddings(anchorset = seu.anchors, reference = seu20_anno_tnk, query = seu80_tnk, new.reduction.name = "pca")

#seu80_tnk <- ProjectUMAP(query = seu80_tnk, query.reduction = "pca", reference = seu20_anno_tnk, reference.reduction = "pca", reduction.model = "umapTNK", reduction.name = "umapTNK", reduction.key = "umapTNK_")

#seu20_tnk_cells <- rownames(na.omit(seu20_anno@reductions$umapTNK))
#celltypes.tmp <- unique(subset(seu20_anno, cells = cells.tmp)$celltype)

# Add UMAP coordinates to metadata
#seu20_anno$UMAP_1 <- seu20_anno@reductions$umap@cell.embeddings[,1]
#seu20_anno$UMAP_2 <- seu20_anno@reductions$umap@cell.embeddings[,2]
#seu80$UMAP_1 <- seu80@reductions$ref.umap@cell.embeddings[,1]
#seu80$UMAP_2 <- seu80@reductions$ref.umap@cell.embeddings[,2]

  
# Now the only difference is the TNK UMAP coordinates. We don't have those for seu80, and perhaps we won't need them
setdiff(colnames(seu20_anno@meta.data), colnames(seu80@meta.data))
setdiff(colnames(seu80@meta.data), colnames(seu20_anno@meta.data))
seu20_anno@reductions
seu80@reductions








# Add complete UMAP coordinates
umap_reduction_df <- data.frame(UMAP_1 = NA, UMAP_2 = NA, rownames = colnames(seu_merge_anno))
umap_reduction_df[colnames(seu20_anno),"UMAP_1"] <- seu20_anno@reductions$umap@cell.embeddings[,1]
umap_reduction_df[colnames(seu20_anno),"UMAP_2"] <- seu20_anno@reductions$umap@cell.embeddings[,2]
umap_reduction_df[colnames(seu80),"UMAP_1"] <- seu80@reductions$ref.umap@cell.embeddings[,1]
umap_reduction_df[colnames(seu80),"UMAP_2"] <- seu80@reductions$ref.umap@cell.embeddings[,2]  
seu_merge_anno[["umap"]] <- CreateDimReducObject(embeddings = umap_reduction_df, key = "umap_")

# Add TNK UMAP coordinates
tnk_reduction_df <- data.frame(umapTNK_1 = NA, umapTNK_2 = NA, rownames = colnames(seu_merge_anno))
tnk_reduction_df[colnames(seu20_T),"umapTNK_1"] <- seu20_T@reductions$umap@cell.embeddings[,1]
tnk_reduction_df[colnames(seu20_T),"umapTNK_2"] <- seu20_T@reductions$umap@cell.embeddings[,2]
tnk_reduction_df[colnames(seu80_TNK),"umapTNK_1"] <- seu80_TNK@reductions$umap@cell.embeddings[,1]
tnk_reduction_df[colnames(seu80_TNK),"umapTNK_2"] <- seu80_TNK@reductions$umap@cell.embeddings[,2]
seu_merge_anno[["umapTNK"]] <- CreateDimReducObject(embeddings = tnk_reduction_df, key = "umapTNK_")

# Check
seu_merge_anno@reductions # should contain umap and umapTNK
apply(seu_merge_anno@reductions$umap, 2, function(x) sum(is.na(x))) # should be 0, 0

nrow(na.omit(seu_merge_anno@reductions$umapTNK)) # should be the same as:
seu_merge_anno@meta.data %>% filter(celltypes %in% tnk_celltypes) %>% nrow







# Load colors from 2.3_PvG-Colors.R
celltype_colors_df <- read.table("../celltype_colors.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE, comment.char = "")
celltype_colors <- setNames(celltype_colors_df$color, celltype_colors_df$celltype)

# Use ggplot to create a UMAP (DimPlot doesn't work b/c the reductions were lost during merge):
seu_merge_anno@meta.data %>%
  sample_frac(1) %>%  # Randomly shuffle rows
  ggplot(aes(x = UMAP_1, y = UMAP_2, color = celltype)) +
  geom_scattermore(pointsize = 8, pixels = c(4096, 4096)) +
  scale_color_manual(values = celltype_colors) +
  theme_bw() +
  theme(aspect.ratio = 1,
        panel.grid = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 3)))

# Save pdf
ggsave("2.5_UMAP_all_cells.pdf", width = 8, height = 4.5)

# Save to persistent disk
saveRDS(seu_merge_anno, "~/250412_SeuratAnno.rds")
#seu_merge <- readRDS("~/250127_seurat_annotated_497K.rds")

# Copy to bucket
#system("gsutil cp ~/250127_seurat_annotated_497K.rds gs://fc-3783b423-62ac-4c69-8c2f-98cb0ee4503b/")

# Compare
seu_old <- readRDS("~/250127_seurat_annotated_497K.rds")

# Basic features
seu_merge
seu_old

# New one has additional "cohort_detail" column
seu_old@meta.data %>% head
seu@meta.data %>% head
setdiff(colnames(seu_old@meta.data), colnames(seu_merge@meta.data))
setdiff(colnames(seu_merge@meta.data), colnames(seu_old@meta.data))

# Cells and cell type are the same
identical(colnames(seu_old), colnames(seu_merge))
identical(seu_old$celltype, seu_merge$celltype)





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
