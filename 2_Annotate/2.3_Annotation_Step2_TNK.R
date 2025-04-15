# Nurefsan Sariipek and Peter van Galen, updated 250410
# Annotate T cells in a more granular way

# Load the libraries
library(tidyverse)
library(Seurat)
library(R.utils)
#library(harmony)
#library(janitor)
#library(ggpubr)

# Start with a clean slate
rm(list=ls())

# Set working directory 
setwd("~/TP53_ImmuneEscape/2_Annotate/")


# 1. DIMENSIONALITY REDUCTION AND CLUSTER --------------------------------------

# Load the Seurat object from 2.1_Subset20percent_and_cluster.R that contains 20% of the cells
seu20 <- readRDS("~/250410_SubsettedSeuratObject.rds")

# Add cell type annotations from 2.2_AnnotationStage1.R
anno_tib <- read_csv("2.2_Step1_celltype_annotations.csv.gz")
anno_df <- data.frame(celltype = anno_tib$celltype, row.names = anno_tib$cell)
seu20 <- AddMetaData(seu20, anno_df)

# Subset to T and NK cells from annotating 20% of the cells
seu20_T <- subset(seu20, celltype == "T cells")

# Cluster and run UMAP. We decided not to rerun VariableFeatures and PCA
seu20_T <- FindNeighbors(seu20_T, dims = 1:20)
seu20_T <- FindClusters(seu20_T, resolution = 1.1)
seu20_T <- RunUMAP(seu20_T, dims = 1:20, return.model = T)

pdf("2.3.1_Step2_UMAP_clusters.pdf", width = 6, height = 6)
DimPlot(seu20_T, group.by = "patient_id", shuffle = T, raster = T) + theme(aspect.ratio = 1, legend.position = "none")
DimPlot(seu20_T, group.by = "seurat_clusters", label = T, raster = T) + theme(aspect.ratio = 1)
FeaturePlot(seu20_T, features = "CD4", reduction = "umap", raster = T) + theme(aspect.ratio = 1)
FeaturePlot(seu20_T, features = "CD8A", reduction = "umap", raster = T) + theme(aspect.ratio = 1)
FeaturePlot(seu20_T, features = "NCAM1", reduction = "umap", raster = T) + theme(aspect.ratio = 1)
dev.off()


# 2. EVALUATE MARKER GENES AND SIGNATURES --------------------------------------

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
# Exhausted (Penter 2021)
FeaturePlot(seu20_T, features = c("PDCD1", "HAVCR2", "LAG3", "TIGIT", "KLRG1", "CD38", "TOX", "TBX21", "EOMES", "CD244", "TNFRSF9", "GZMB", "PRF1"), shuffle = T)

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


# 3. # ANNOTATE AND SAVE -------------------------------------------------------

# Add T cell cluster names
Tcell.cluster.ids <- c("CD4 Memory","CD8 Memory","CD4 Naive","Treg","CD8 Effector","CD4 Effector Memory","Gamma-Delta T","CD8 Exhausted","CD56 Dim NK","CD8 Naive","NK-T","CD56 Dim NK","CD56 Bright NK", "Adaptive NK")
names(Tcell.cluster.ids) <- levels(seu20_T)
seu20_T <- RenameIdents(seu20_T, Tcell.cluster.ids)
seu20_T$celltype <- Idents(seu20_T)

DimPlot(seu20_T, label = T) + theme(aspect.ratio = 1, legend.position = "none")

# Save plot
p1 <- DimPlot(seu20_T, reduction = "umap", repel = T, label = T, raster = T) +
  theme(aspect.ratio = 1, legend.position = "none")

pdf("2.3.2_Step2_Annotated_clusters.pdf", width = 6, height = 6)
print(p1)
dev.off()

# Save a tibble with all cell type annotations and predictions
celltypes_tib <- as_tibble(seu20_T@meta.data, rownames = "cell") %>% select(cell, celltype)
write_csv(celltypes_tib, file = "2.3_Step2_TNK_annotations.csv")
gzip("2.3_Step2_TNK_annotations.csv", overwrite = T, remove = T)

# Save a tibble with all UMAP coordinates and projections
coordinates_tib <- as_tibble(seu20_T@reductions$umap@cell.embeddings, rownames = "cell")
write_csv(coordinates_tib, file = "2.3_Step2_umapTNK_coordinates.csv")
gzip("2.3_Step2_umapTNK_coordinates.csv", overwrite = T, remove = T)


