# Nurefsan Sariipek and Peter van Galen, updated 250410
# Annotate T cells in a more granular way

# Load the libraries
library(tidyverse)
library(Seurat)
library(harmony) # only needed for section 1
#library(janitor)
#library(ggpubr)

# Start with a clean slate
rm(list=ls())

# Set working directory 
setwd("~/TP53_ImmuneEscape/2_Annotate/")


# 1. DIMENSIONALITY REDUCTION AND CLUSTER --------------------------------------

# Load the data from 1.1_CreateSeuratObject.R
seu <- readRDS(file = "~/250409_MergedSeuratObject.rds")

# Load the Seurat object from 2.1_Subset20percent_and_cluster.R that contains 20% of the cells (only needed for subsetting)
seu20 <- readRDS("~/250410_SubsettedSeuratObject.rds")

# Add cell type annotations from 2.2_AnnotationStage1.R
anno_tib <- read_csv("2.2_Stage1_celltype_annotations.csv.gz")
anno_df <- data.frame(celltype = anno_tib$celltype, row.names = anno_tib$cell)
seu <- AddMetaData(seu, anno_df)

# Subset to T and NK cells from annotating 20% of the cells
seu20_fresh <- subset(seu, cells = colnames(seu20))
seu20_T <- subset(seu20_fresh, celltype == "T cells")

# Then perform standard clustering and UMAP
seu20_T <- NormalizeData(seu20_T)
seu20_T <- FindVariableFeatures(seu20_T)
p1 <- VariableFeaturePlot(seu20_T)
LabelPoints(plot = p1, points = head(VariableFeatures(seu20_T), 20), xnudge = 0, ynudge = 0) +
  theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5)) +
  ggtitle("Variable genes")

# Dimensionality reduction
seu20_T <- ScaleData(seu20_T, features = rownames(seu20_T))
seu20_T <- RunPCA(seu20_T, features = VariableFeatures(object = seu20_T))
seu20_T <- RunHarmony(object = seu20_T, group.by.vars = c("patient_id"), plot_convergence = T)
ElbowPlot(seu20_T, reduction = "harmony")
seu20_T <- RunUMAP(seu20_T, reduction = "harmony", dims = 1:18, return.model = T)

# Cluster data
seu20_T <- FindNeighbors(seu20_T, reduction = "harmony", dims = 1:18)
seu20_T <- FindClusters(seu20_T, resolution = 0.5)

# Visualize different UMAPs
DimPlot(seu20_T, reduction = "umap", group.by = "patient_id", shuffle = T) +
  theme(aspect.ratio = 1, legend.position = "none")
DimPlot(seu20_T, reduction = "umap", group.by = "seurat_clusters", label = T) + theme(aspect.ratio = 1)
FeaturePlot(seu20_T, features = "CD8A", reduction = "umap") + theme(aspect.ratio = 1)
FeaturePlot(seu20_T, features = "CD4", reduction = "umap") + theme(aspect.ratio = 1)
FeaturePlot(seu20_T, features = "NCAM1", reduction = "umap") + theme(aspect.ratio = 1)

# Save for faster loading next time
seu20_T@assays$RNA@layers$scale.data <- NULL
saveRDS(seu20_T, "~/250413_SubsettedSeuratTNK_clusters.rds")


# 2. VISUALIZE MARKER GENES AND SIGNATURES -------------------------------------

# Load object with seurat_clusters
seu20_T <- readRDS("~/250413_SubsettedSeuratTNK_clusters.rds")

# Look at clusters
DimPlot(seu20_T, reduction = "umap", group.by = "seurat_clusters", label = T) + theme(aspect.ratio = 1)

# Run Find Markers, convert to a tibble and save it as a csv
Tcell_markers <- FindAllMarkers(seu20_T, min.pct = .3, logfc.threshold = .3)
top50_markers_tib <- as_tibble(Tcell_markers) %>% filter(avg_log2FC > 0) %>%
  group_by(cluster) %>% arrange(-avg_log2FC) %>% slice_head(n = 50) %>%
  select(cluster, gene) %>%
  pivot_wider(names_from = "cluster", values_from = "gene",
              names_prefix = "cluster_", values_fn = list) %>%
  unnest(cols = everything())
write_csv(top50_markers_tib, file = "2.3_Tcell_markers.csv")

# Evaluate T cell features to distinguish populations
FeaturePlot(seu20_T, features = c("CD8A", "CD8B", "CD4", "NCAM1", "IL10", "TGFB",
                                "GATA3", "TCF7", "SELL", "CCR7", "SELL", "TMIGD2",
                                "LEF1", "CD28", "CD27"))
FeaturePlot(seu20_T, features = c("CD8A", "CD8B", "CD4", "NCAM1", "IFNG", "CCL3",
                                "PDCD1", "TIGIT", "LAG3", "HAVCR2", "CTLA4"))

# Naive
FeaturePlot(seu20_T, features = c("CCR7", "TCF7", "LEF1", "SELL"))

# CD8 Naive 
FeaturePlot(seu20_T, features = c("CD8A", "CD8B", "LEF1", "SELL", "CCR7",
                                "TCF7", "CD27", "CD28", "S1PR1"))

# Effector 
FeaturePlot(seu20_T, features = c("GZMA", "GZMB", "PRF1", "CX3CR1"))

# Cytotoxicity 
FeaturePlot(seu20_T, features = c("CX3CR1", "GMZH", "GNLY", "FGFBP2"))

# Cytotoxic
FeaturePlot(seu20_T, features = c("FCGR3A", "KLRG1", "PRF1", "GZMB"))

# TIL atlas
# CD8 tpex
FeaturePlot(seu20_T, features =c("LAG3", "XCL1", "CRTAM", "IFNG", "CCL4", "PDCD1", "DUSP4", "CD8A", "ZEB2", "NR4A2", "SLA", "NKG7", "TIGIT", "CTSW", "TNFRSF9", "TOX", "LYST", "TNFSF4", "CCL3", "GZMB", "RAB27A", "PRF1", "CD70", "PLSCR1", "CXCL13"))

# CD4 naive like
FeaturePlot(seu20_T, features = c("CCR7", "SELL", "CD40LG", "IL7R", "TCF7", "LEF1", "GPR183", "KLRB1", "LTB", "MAL", "PASK", "AQP3", "TRAT1"))

# Cytotoxic
FeaturePlot(seu20_T, features = c("NKG7", "CCL4", "CST7", "GZMA", "GZMB", "IFNG", "CCL3"))

# Effector memory
FeaturePlot(seu20_T, features = c("KLRG1", "LYAR", "GZMK", "GZMM", "TXNIP", "FCRL6"))

# Memory 
FeaturePlot(seu20_T, features = c("TCF7", "CCR7", "IL7R", "SELL"))

# T-NK
FeaturePlot(seu20_T, features = c("TCF7", "GZMM", "KLRC3", "KLRB1", "GNLY"))

# Exhausted
FeaturePlot(seu20_T, features = c("LAYN", "ITGAE", "PDCD1", "CTLA4", "HAVCR2", "LAG3", "TIGIT"))

# Terminally Exhausted
FeaturePlot(seu20_T, features = c("TRAV24", "TRBV5-5", "CCDC50", "HLA-DRB1", "VCAM1", "CTSW", "NKG7", "CXCL13", "HLA-DPA1", "DUSP4", "GZMB", "HLA-DRA"))

# CD8 T Effector Memory cells
FeaturePlot(seu20_T, features = c("S1PR1", "GPR183", "CCR7", "ANXA1", "TCF7", "IL7R", "MBP", "VIM", "GYG1", "GPR183"))

# CD8 T progenitors
FeaturePlot(seu20_T, features = c("CAV1", "GNG4", "XCL1", "CRTAM", "CXCL13", "GEM", "XCL2"))

# TREG
FeaturePlot(seu20_T, features = c("FOXP3", "IL2RA", "SELL", "TNFRSF4"))

# MAIT cells
FeaturePlot(seu20_T, features = c("KLRB1", "NCR3", "ZBTB16", "SLC4A10", "RORC"))

# Heatmap the genes in your interest, providing them as gene_list
# For the genes down below, we have used various published datasets, which will be referenced in our paper.

# All genes (from Melanoma paper)
gene_list <- c("CD3E", "CD4", "CD8A", "SELL", "CCR7", "IL7R", "CD28", "FAS", "CD27", "ITGAE", "ITGAL", "ITGAM", "ITGAX", "PDCD1", "TIGIT", "HAVCR2", "LAG3", "CTLA4", "VTCN1", "CD244", "KLRG1", "TNFRSF14", "BTLA", "CD160", "CD38", "ENTPD1", "NT5E", "CD69", "IL2RA", "ICOS", "TNFRSF4", "TNFRSF9", "HLA-DRA", "CD40LG", "GZMA", "GZMB", "GZMH", "GZMK", "GZMM", "PRF1", "NKG7", "GNLY", "IFNG", "FASLG", "TNF", "IL17A", "IL2", "LEF1", "TCF7", "EOMES", "TBX21", "PRDM1", "TOX", "GATA3", "ID2", "ID3", "NR4A1", "ZNF683", "FOXP3", "MKI67", "TOP2A", "TRGV9", "TRDV2", "KLRB1", "KLRC3")

# CD8 genes
gene_list <- c("CD8", "CD8A", "CD8B", "CD4", "SELL", "ILR7A", "CD28", "CD27", "TIGIT", "HAVCR2", "LAG3", "CTLA4", "KIRG2" , "PDCD1", "CD160", "CD38", "ENTPD1", "IL2RA", "ICOS", "HLA-DRA", "CD40LG", "NKG7", "GNLY", "CST7", "PRF1", "GZMK", "GZMH", "GZMA", "GZMB", "IFNG", "TNF", "IL17A", "NKG7", "GNLY", "FASLG", "TRGV9", "TRDV2", "KLRB1", "KLRC3")

# CD4 genes
gene_list <- c("CD4", "NKG7", "GNLY", "CCL4", "CST7", "PRF1", "GZMK", "GZMH", "GZMA", "GZMB", "IFNG", "CCL3", "PDCD1", "TIGIT", "LAG3", "HAVCR2", "CTLA4", "FOXP3", "IKZF2", "TIMP1", "CCL5", "ANXA1", "S100A11", "ANXA2", "IL2RA", "CCR10", "TNFRSF4", "TIMP4", "CCR7", "SELL", "TCF7", "LEF1", "CD27", "CD28")

# CD56 NK cells genes
gene_list <- c("NCAM", "CD56", "CD57", "CD16", "FCGR3A", "KIR", "NKG2C", "NKG2A", "CD57", "GZMK", "XCL1", "ITGA1", "IL7R", "TCF7", "GPR183", "GZMB", "PRF1", "CX3CR1", "CD7", "FCER1G", "RUNX2", "BACH2", "MAHML3", "TCF4", "ZEB1", "ZMAT4", "ZBTB16", "ZNF516", "ZEB2", "ARID5B", "KLRB1", "KLRC2", "CD3E", "PATL2", "ZBTB38")

# AverageExpression will normalize and scale the data when you set the return.seurat = T
seu20_T_avg <- AggregateExpression(seu20_T, return.seurat = T)
seu20_T_avg_cd8 <- AggregateExpression(seu20_T %>% subset(seurat_clusters %in% c(1,4,6,7)), return.seurat = T)
seu20_T_avg_nk <- AggregateExpression(seu20_T %>% subset(seurat_clusters %in% c(8,10,11,12,13)), return.seurat = T)
seu20_T_avg_cd4 <- AggregateExpression(seu20_T %>% subset(seurat_clusters %in% c(0,2,3)), return.seurat = T)

# Visualize it as a heatmap with the gene_list you provide, change the object to visualize different ones
DoHeatmap(seu20_T_avg_nk, features = gene_list, draw.lines = FALSE)  +  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 5, name = "RdBu"))) + theme_pubr(base_size = 10) + guides(colour = "none")

# Also visualize signature module scores from others, starting with Kyle Romine's
kyle_markerGenes <- read.csv(file = "Signatures/KR_TCellTypeSignatures.csv", header = T)

# Add module scores 
for (n in names(kyle_markerGenes)) {
  print(n)
  seu20_T <- AddModuleScore(object = seu20_T, features = kyle_markerGenes[n], name = n)
  colnames(seu20_T@meta.data) <- gsub(str_c(n, "1$"), str_c(n, "_Score"), colnames(seu20_T@meta.data))
  }

# Visualize
FeaturePlot(seu20_T, features = c("CD4_Naïve_Score", "CD56_dim_NK_Score", "CD8_Term_Eff_Score", "CD8_GZMK_Exh_Score", "CD8_EM_Score", "CD8_Naïve_Score", "NKT_Score", "CD4_CM_Score", "MAIT_Score", "Tregs_Score", "CD56_Bright_NK_Score"), raster = T)
ggsave("FeaturePlots/Tsubset_Kyle_FeaturePlots.pdf", width = 10, height = 8)


# 3. ADD ANNOTATIONS -----------------------------------------------------------

# Based on the visualizations generated above, add T cell cluster names
Tcell.cluster.ids <- c("CD4 Memory", "CD8 Memory", "CD4 Naïve", "Treg", "CD8 Effector", "CD4 Effector Memory", "γδ T", "CD8 Exhausted", "CD56 Dim NK", "CD8 Naïve", "NK T", "CD56 Dim NK", "CD56 Bright NK", "Adaptive NK")
names(Tcell.cluster.ids) <- levels(seu20_T)
seu20_T <- RenameIdents(seu20_T, Tcell.cluster.ids)
seu20_T$celltype <- Idents(seu20_T)
# See the levels
levels(seu20_T$celltype)

# See the new annotated UMAP
DimPlot(seu20_T, reduction = "umap", repel = T, group.by = "celltype", label = T) + theme(aspect.ratio = 1)

# Save for projection of remaining T cells in 2.3_Complete_SeuratAnno.R (possibly obsolete)
seu20_T@assays$RNA@layers$scale.data <- NULL
seu20_T$RNA_snn_res.1 <- NULL
seu20_T$RNA_snn_res.1.1 <- NULL
seu20_T$seurat_clusters <- NULL
seu20_T@meta.data[,grepl("_Score", colnames(seu20_T@meta.data))] <- NULL
saveRDS(seu20_T, file = "250411_AnnotatedSeurat_TNK.rds")

# Add merged cell type annotations to Seurat object
meta <- seu20@meta.data
tnk_meta <- seu20_T@meta.data
meta$celltype <- as.character(meta$celltype)
tnk_meta$celltype <- as.character(tnk_meta$celltype)
meta$cell <- rownames(meta)
tnk_meta$cell <- rownames(tnk_meta)
meta$celltype[meta$cell %in% tnk_meta$cell] <- tnk_meta$celltype
seu20@meta.data <- meta

# Turn celltype into a factor and add levels 
seu20@meta.data$celltype <- factor(seu20@meta.data$celltype,levels = c("Progenitors", "Early Erythroids", "Mid Erythroids", "Late Erythroids", "Pro Monocytes", "Monocytes", "Non Classical Monocytes", "cDC",  "pDC", "Pro B cells", "Pre-B", "B cells", "Plasma cells", "CD4 Naïve", "CD4 Effector Memory", "CD4 Memory", "Treg", "CD8 Naïve", "CD8 Effector", "CD8 Memory", "CD8 Exhausted",  "γδ T", "NK T", "Adaptive NK", "CD56 Bright NK", "CD56 Dim NK", "Cycling T-NK cells", "UD1", "UD2", "UD3"))
seu20@active.ident <- seu20$celltype 
DimPlot(seu20, reduction = "umap", repel = T, label = T) + theme(aspect.ratio = 1)

# Add T-NK umap coordinates to the merged object (likely obsolete)
umapTNK_df <- seu20_T@reductions$umap@cell.embeddings
colnames(umapTNK_df) <- paste0("umapTNK_", 1:2)
seu20[["umapTNK"]] <- CreateDimReducObject(embeddings = umapTNK_df, key = "umapTNK_")

# Remove unnecessary information
seu20@assays$RNA@layers$scale.data <- NULL
seu20$RNA_snn_res.1 <- NULL
seu20$seurat_clusters <- NULL
seu20$cell <- NULL

# Save the merged seu as a new object
saveRDS(seu20, file = "~/250411_SubsettedAnnotatedSeuratObject.rds")



# Compare new Seurat object with a previous iteration
seu_old <- readRDS("~/250127_seu_annotated_merged_no-scale.rds")

# Basic features - same
seu_old
seu20

# Some differences in metadata
setdiff(colnames(seu_old@meta.data), colnames(seu20@meta.data))
setdiff(colnames(seu20@meta.data), colnames(seu_old@meta.data))

# Reductions
seu_old@reductions
seu20@reductions

# Cells and cell type are the same
identical(colnames(seu_old), colnames(seu20))
identical(seu_old$celltype, seu20$celltype)




# PROJECT umapTNK FOR ALL OTHER CELLS ------------------------------------------

tnk_celltypes <- c("CD4 Naïve", "CD4 Memory", "CD4 Effector Memory", "Treg", "CD8 Naïve", "CD8 Effector", 
                   "CD8 Memory", "CD8 Exhausted", "γδ T", "NK T", "Adaptive NK", "CD56 Bright NK", "CD56 Dim NK")

# WIP: run the first section under ANNOTATE AND SAVE
# TO TEST: make 10x smaller
#seu20_T_backup <- seu20_T
#seu20_T <- seu20_T_backup
#seu20_T <- subset(seu20_T, cells = colnames(seu20_T)[seq(1, ncol(seu20_T), by = 10)])
#seu20_T$celltype <- factor(seu20_T$celltype, levels = tnk_celltypes)

seu_merge <- readRDS(file = "~/250409_MergedSeuratObject.rds")
predictions_tib <- read_csv("2.3_CellTypePredictions.csv")
predictions_df <- data.frame(celltype = predictions_tib$celltype, row.names = predictions_tib$cell)
seu_merge <- AddMetaData(seu_merge, predictions_df)
seu80_T <- subset(seu_merge, celltype %in% tnk_celltypes)

# UMAP projection of the remaining 80% of T/NK cells (similar to https://satijalab.org/seurat/articles/integration_mapping.html#unimodal-umap-projection)
seu20_T@assays$RNA@layers$scale.data <- NULL
seu80_T <- NormalizeData(seu80_T)
seu.anchors <- FindTransferAnchors(reference = seu20_T, query = seu80_T, dims = 1:20, reference.reduction = "pca")
seu80_T <- IntegrateEmbeddings(anchorset = seu.anchors, reference = seu20_T, query = seu80_T, new.reduction.name = "ref.pca")
seu80_T <- ProjectUMAP(query = seu80_T, query.reduction = "ref.pca", reference = seu20_T, reference.reduction = "pca", reduction.model = "umap", reduction.name = "ref.umap", reduction.key = "refUMAP_")

DimPlot(seu20_T, group.by = "patient_id", reduction = "umap") + theme(aspect.ratio = 1)
DimPlot(seu80_T, group.by = "patient_id", reduction = "ref.umap") + theme(aspect.ratio = 1)

# I've now tried a number of strategies to no avail.


