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
seu <- readRDS(file = "~/250409_MergedSeuratObject.rds")
seu20 <- readRDS("~/250410_SubsettedSeuratObject.rds")

# Remove the annotated cells from the full object
seu80 <- subset(seu, cells = setdiff(colnames(seu), colnames(seu20)))
# Check that the cells numbers make sense
identical(ncol(seu80), ncol(seu)-ncol(seu20))

# Add cell type annotations: start with all cells
anno1_tib <- read_csv("2.2_Step1_celltype_annotations.csv.gz")
anno1_df <- data.frame(celltype = anno1_tib$celltype, row.names = anno1_tib$cell)
# TNK cells
anno2_tib <- read_csv("2.3_Step2_TNK_annotations.csv.gz")
anno2_df <- data.frame(celltype = anno2_tib$celltype, row.names = anno2_tib$cell)
# Replace "T cells" with granular annotations
anno1_df[rownames(anno1_df) %in% rownames(anno2_df),"celltype"] <- anno2_df$celltype
# Add to Seurat object
seu20 <- AddMetaData(seu20, anno1_df)


# PROJECT FULL UMAP AND PREDICT CELL TYPES -------------------------------------

# UMAP projection of the remaining 80% of cells (similar to https://satijalab.org/seurat/articles/integration_mapping.html#unimodal-umap-projection)
seu80 <- NormalizeData(seu80)
seu.anchors <- FindTransferAnchors(reference = seu20, query = seu80, dims = 1:18, reference.reduction = "pca")
seu80 <- IntegrateEmbeddings(anchorset = seu.anchors, reference = seu20, query = seu80, new.reduction.name = "ref.pca")
seu80 <- ProjectUMAP(query = seu80, query.reduction = "ref.pca", reference = seu20, reference.reduction = "pca", reduction.model = "umap")

# Cell type prediction for the remaining 80% of cells
predictions_df <- TransferData(anchorset = seu.anchors, refdata = seu20$celltype, dims = 1:18)

# Compare original and projected/predicted data
colnames(predictions_df) <- gsub("predicted.id", "celltype", colnames(predictions_df))

# Compare UMAPs
seu80 <- AddMetaData(seu80, metadata = select(predictions_df, celltype))
p1 <- DimPlot(seu20, reduction = "umap", group.by = "celltype", shuffle = T, label = T, raster = T) +
  theme(aspect.ratio = 1, legend.position = "none")
p2 <- DimPlot(seu80, reduction = "ref.umap", group.by = "celltype", shuffle = T, label = T, raster = T) +
  theme(aspect.ratio = 1, legend.position = "none")

pdf("2.4.1_Projections_and_predictions.pdf", width = 15, height = 10)
p1 + p2
dev.off()


# MERGE AND SAVE ---------------------------------------------------------------

# Add complete UMAP coordinates
umap_reduction_df <- data.frame(UMAP_1 = rep(NA, ncol(seu)), UMAP_2 = rep(NA, ncol(seu)), row.names = colnames(seu))
umap_reduction_df[colnames(seu20),"UMAP_1"] <- seu20@reductions$umap@cell.embeddings[,1]
umap_reduction_df[colnames(seu20),"UMAP_2"] <- seu20@reductions$umap@cell.embeddings[,2]
umap_reduction_df[colnames(seu80),"UMAP_1"] <- seu80@reductions$ref.umap@cell.embeddings[,1]
umap_reduction_df[colnames(seu80),"UMAP_2"] <- seu80@reductions$ref.umap@cell.embeddings[,2]  
seu[["umap"]] <- CreateDimReducObject(embeddings = as.matrix(umap_reduction_df), key = "umap_")
DimPlot(seu, reduction = "umap", group.by = "patient_id") + theme(aspect.ratio = 1)

# Add TNK UMAP coordinates. The projection for additional T cells does not work
coord_tib <- read_csv("2.3_Step2_umapTNK_coordinates.csv.gz")
coord_mat <- as.matrix(data.frame(coord_tib[,"umap_1"], coord_tib[,"umap_2"], row.names = coord_tib$cell))
seu[["umapTNK"]] <- CreateDimReducObject(embeddings = coord_mat, key = "umapTNK_")
DimPlot(seu, reduction = "umapTNK", group.by = "patient_id") + theme(aspect.ratio = 1)

# Add cell type annotations
celltypes_df <- rbind(seu20@meta.data[,"celltype",drop = F],
                      seu80@meta.data[,"celltype",drop = F])
seu <- AddMetaData(seu, celltypes_df)

# Set cell types as a logically ordered factor
celltypes <- c("Progenitors", "Early Erythroid", "Mid Erythroid", "Late Erythroid", "Pro Monocytes",
               "Monocytes", "Non-Classical Monocytes", "cDC",  "pDC", "Pro-B", "Pre-B", "B cells",
               "Plasma cells", "CD4 Naive", "CD4 Memory", "CD4 Effector Memory", "Treg",
               "CD8 Naive", "CD8 Memory", "CD8 Effector", "CD8 Exhausted", "Gamma-Delta T", "NK-T",
               "Adaptive NK", "CD56 Bright NK", "CD56 Dim NK", "Cycling T-NK",
               "UD1", "UD2", "UD3")
# Check that the overlap is perfect
all(seu$celltype %in% celltypes)
all(celltypes %in% seu$celltype)
seu$celltype <- factor(seu$celltype, levels = celltypes)
seu@active.ident <- seu$celltype

# Add normalized data slot
seu <- NormalizeData(seu)

# Some checks
seu@reductions # should contain umap and umapTNK
apply(seu@reductions$umap@cell.embeddings, 2, function(x) sum(is.na(x))) # should be 0, 0
nrow(na.omit(seu@reductions$umapTNK@cell.embeddings)) # should be 49745
tnk_celltypes <- c("CD4 Naive", "CD4 Memory", "CD4 Effector Memory", "Treg",
                   "CD8 Naive", "CD8 Memory", "CD8 Effector", "CD8 Exhausted", "Gamma-Delta T", "NK-T",
                   "Adaptive NK", "CD56 Bright NK", "CD56 Dim NK")
seu@meta.data %>% filter(celltype %in% tnk_celltypes) %>% nrow # should be 248936

# Load colors from 2_Annotate/Colors.R
celltype_colors_df <- read.table("../celltype_colors.txt", sep = "\t", header = T, stringsAsFactors = F, comment.char = "")
celltype_colors <- setNames(celltype_colors_df$color, celltype_colors_df$celltype)

# Show UMAP for all cells
DimPlot(seu, reduction = "umap", shuffle = T, raster = T, cols = celltype_colors) +
  theme_bw() +
  theme(aspect.ratio = 1, panel.grid = element_blank())

ggsave("2.4.1_All_cells_annotated.pdf", width = 8, height = 4.5)

# Show UMAP for TNK cells. Only 20% of the cells have UMAP coordinates.
DimPlot(seu, reduction = "umapTNK", shuffle = T, raster = T, pt.size = 2, cols = celltype_colors) +
  theme_bw() +
  theme(aspect.ratio = 1, panel.grid = element_blank())

ggsave("2.4.2_TNK_cells_annotated.pdf", width = 6, height = 4)

# Save to persistent disk
saveRDS(seu, "~/250415_Seurat_all_cells_annotated.rds")

# Copy to bucket (Terminal)
#system("gsutil cp ~/250415_Seurat_all_cells_annotated.rds gs://fc-3783b423-62ac-4c69-8c2f-98cb0ee4503b/")



# Compare with previous clustering
seu_old <- readRDS("~/250127_seurat_annotated_497K.rds")

# Basic features
seu
seu_old

# New one has additional "cohort_detail" column
seu_old@meta.data %>% head
seu@meta.data %>% head
setdiff(colnames(seu_old@meta.data), colnames(seu@meta.data))
setdiff(colnames(seu@meta.data), colnames(seu_old@meta.data))

# Check cells and cell types
identical(colnames(seu_old), colnames(seu))
identical(sort(colnames(seu_old)), sort(colnames(seu)))
celltypes_tib <- select(as_tibble(seu_old@meta.data, rownames = "cell"), cell, celltype) %>%
  full_join(., select(as_tibble(seu@meta.data, rownames = "cell"), cell, celltype), by = "cell") %>%
  group_by(celltype.x, celltype.y) %>% count
view(celltypes_tib)



