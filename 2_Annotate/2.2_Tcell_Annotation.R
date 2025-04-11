# Nurefsan Sariipek and Peter van Galen, updated 250410
# Annotate T cells in a more granular way

# Load the libraries
library(tidyverse)
library(Seurat)
library(janitor)
library(ggpubr)

# Start with a clean slate
rm(list=ls())

# Set working directory 
setwd("~/TP53_ImmuneEscape/2_Annotate/")

# Load the saved Seurat object from the last step that contains only 20% of the cells
seu <- readRDS("~/250410_SubsettedSeuratObject.rds")

# Add cell type annotations as Seurat identity and metadata column as indicated at the end of 2.1_Annotate.R
Idents(seu) <- "seurat_clusters"
seu.cluster.ids <- c("T cells", "T cells", "Monocytes", "T cells", "Late Erythroids", "T cells", "Mid Erythroids", "UD1", "B cells", "Monocytes", "Non Classical Monocytes", "T cells", "Pro Monocytes", "Late Erythroids", "Early Erythroids", "UD3", "Early Erythroids", "Pre-B","Progenitors", "cDC", "Pro B cells", "UD1", "Early Erythroids", "Plasma cells", "Cycling T-NK cells", "pDC", "Progenitors", "Progenitors", "Pro B cells", "Late Erythroids", "UD2", "UD2", "UD2")
names(seu.cluster.ids) <- levels(seu)
seu <- RenameIdents(seu, seu.cluster.ids)
seu@meta.data$celltype = Idents(seu)

# Subset only T cells from all metadata
seu_T <- subset(x = seu, subset = seurat_clusters %in% c(0,1,3,5,11))

# Run FindVariables again (replace the ones inherited from the complete object)
seu_T <- FindVariableFeatures(seu_T)

# Decide dimensions according to elbow plot below
ElbowPlot(seu_T)
seu_T <- FindNeighbors(seu_T, dims = 1:20)
seu_T <- FindClusters(seu_T, resolution = 1.1)

# Run UMAP
seu_T <- RunUMAP(seu_T, dims = 1:20)

# Visualize the UMAP
umap_plot <- DimPlot(seu_T, reduction = "umap", label = TRUE) + theme(aspect.ratio = 1)

# Save the UMAP
pdf("2.2.1_UMAP_Tcells.pdf", width = 10, height = 8)
print(umap_plot)
dev.off()

# Run Find Markers, convert to a tibble and save it as a csv
Tcell_markers <- FindAllMarkers(seu_T, min.pct = .3, logfc.threshold = .3)
Tcell_markers_tibble <- as_tibble(Tcell_markers)
write.csv(Tcell_markers_tibble, file = "2.2_Tcell_markers.csv")

# Assuming the file has columns like 'Subset' and 'Score' to determine top markers,
# modify the column names based on your data.
# For example:
# Subset: The column that identifies different T cell subsets.
# Score: The column based on which you select the top 10 (e.g., p-value, log2FoldChange, etc.).

# Get the top 10 markers per cluster based on fold change
top_markers <- Tcell_markers_tibble %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) %>%
  arrange(cluster, desc(avg_log2FC))

# Save the top markers to a new CSV
output_file <- "2.2_Top10_Tcell_markers.csv"
write.csv(top_markers, output_file, row.names = FALSE)


######### Annotation ##########

# T cell Features to distinguish populations
FeaturePlot(seu_T, features = c("CD8A", "CD8B", "CD4", "NCAM1","IL10","TGFB","GATA3","TCF7","SELL","CCR7","SELL","TMIGD2","LEF1","CD28","CD27"))
FeaturePlot(seu_T, features = c("CD8A", "CD8B", "CD4", "NCAM1","IFNG", "CCL3", "PDCD1", "TIGIT", "LAG3", "HAVCR2", "CTLA4"))

# Naive
FeaturePlot(seu_T, features = c("CCR7","TCF7","LEF1", "SELL"))

# CD8 Naive 
FeaturePlot(seu_T, features = c("CD8A", "CD8B", "LEF1", "SELL", "CCR7", "TCF7", "CD27", "CD28","S1PR1"))

# Effector 
FeaturePlot(seu_T, features = c("GZMA","GZMB","PRF1","CX3CR1"))

# Cytotoxicity 
FeaturePlot(seu_T, features = c("CX3CR1","GMZH","GNLY","FGFBP2"))

# Cytotoxic
FeaturePlot(seu_T, features = c("FCGR3A", "KLRG1", "PRF1", "GZMB"))

# TIL Atlas
# CD8 tpex
FeaturePlot(seu_T, features =c("LAG3", "XCL1", "CRTAM","IFNG", "CCL4", "PDCD1", "DUSP4", "CD8A", "ZEB2", "NR4A2", "SLA", "NKG7", "TIGIT", "CTSW", "TNFRSF9", "TOX", "LYST", "TNFSF4", "CCL3", "GZMB", "RAB27A", "PRF1", "CD70", "PLSCR1","CXCL13"))

# CD4 naive like
FeaturePlot(seu_T, features = c("CCR7", "SELL", "CD40LG","IL7R", "TCF7", "LEF1", "GPR183", "KLRB1", "LTB", "MAL", "PASK", "AQP3", "TRAT1"))      

# Cytotoxic
FeaturePlot(seu_T, features = c("NKG7", "CCL4", "CST7", "GZMA", "GZMB", "IFNG", "CCL3"))

# Effector memory
FeaturePlot(seu_T, features = c("KLRG1", "LYAR", "GZMK", "GZMM", "TXNIP","FCRL6"))

# Memory 
FeaturePlot(seu_T, features = c("TCF7","CCR7","IL7R","SELL"))

# T-NK
FeaturePlot(seu_T, features = c("TCF7","GZMM","KLRC3","KLRB1","GNLY"))

# Exhausted
FeaturePlot(seu_T, features = c("LAYN", "ITGAE", "PDCD1", "CTLA4", "HAVCR2", "LAG3", "TIGIT"))

# Terminally Exhausted
FeaturePlot(seu_T, features = c("TRAV24","TRBV5-5","CCDC50","HLA-DRB1", "VCAM1", "CTSW", "NKG7", "CXCL13", "HLA-DPA1", "DUSP4", "GZMB","HLA-DRA"))

# CD8 T Effector Memory cells
FeaturePlot(seu_T, features = c("S1PR1", "GPR183", "CCR7", "ANXA1", "TCF7", "IL7R", "MBP", "VIM","GYG1","GPR183"))

# CD8 T progenitors
FeaturePlot(seu_T, features = c("CAV1", "GNG4", "XCL1", "CRTAM", "CXCL13", "GEM","XCL2"))

# TREG
FeaturePlot(seu_T, features = c("FOXP3","IL2RA","SELL","TNFRSF4"))

# MAIT cells
FeaturePlot(seu_T, features = c("KLRB1", "NCR3", "ZBTB16", "SLC4A10", "RORC"))

# Heatmap the genes in your interest, providing them as gene_list
# For the genes down below, we have used various published datasets, which will be referenced in our paper.

# All genes (from Melanoma paper)
gene_list <- c("CD3E","CD4","CD8A","SELL","CCR7","IL7R","CD28","FAS","CD27","ITGAE","ITGAL","ITGAM","ITGAX","PDCD1","TIGIT","HAVCR2","LAG3","CTLA4","VTCN1","CD244","KLRG1","TNFRSF14","BTLA","CD160","CD38","ENTPD1","NT5E","CD69", "IL2RA","ICOS","TNFRSF4","TNFRSF9","HLA-DRA","CD40LG","GZMA","GZMB","GZMH","GZMK","GZMM","PRF1","NKG7","GNLY","IFNG","FASLG","TNF","IL17A","IL2","LEF1","TCF7","EOMES","TBX21","PRDM1","TOX","GATA3","ID2","ID3","NR4A1","ZNF683","FOXP3","MKI67","TOP2A","TRGV9","TRDV2","KLRB1","KLRC3")

# CD8 genes
gene_list <- c("CD8", "CD8A", "CD8B", "CD4", "SELL", "ILR7A", "CD28", "CD27","TIGIT","HAVCR2", "LAG3", "CTLA4", "KIRG2" ,"PDCD1","CD160","CD38","ENTPD1","IL2RA","ICOS","HLA-DRA","CD40LG", "NKG7", "GNLY", "CST7", "PRF1", "GZMK", "GZMH", "GZMA", "GZMB", "IFNG", "TNF", "IL17A", "NKG7", "GNLY", "FASLG", "TRGV9", "TRDV2", "KLRB1", "KLRC3")

# CD4 genes
gene_list <- c("CD4","NKG7","GNLY","CCL4","CST7","PRF1","GZMK","GZMH","GZMA","GZMB","IFNG","CCL3","PDCD1","TIGIT","LAG3","HAVCR2","CTLA4","FOXP3","IKZF2","TIMP1","CCL5","ANXA1","S100A11","ANXA2", "IL2RA","CCR10","TNFRSF4","TIMP4","CCR7","SELL","TCF7","LEF1","CD27","CD28")

# CD56 NK cells genes
gene_list <- c("NCAM","CD56","CD57", "CD16","FCGR3A", "KIR", "NKG2C", "NKG2A", "CD57","GZMK","XCL1", "ITGA1", "IL7R", "TCF7","GPR183","GZMB","PRF1","CX3CR1","CD7","FCER1G", "RUNX2","BACH2","MAHML3","TCF4","ZEB1", "ZMAT4","ZBTB16","ZNF516","ZEB2","ARID5B","KLRB1","KLRC2","CD3E","PATL2", "ZBTB38")

# AverageExpression will normalize and scale the data when you set the return.seurat = T
seu_T_avg <- AggregateExpression(seu_T, return.seurat = T)
seu_T_avg_cd8 <- AggregateExpression(seu_T %>% subset(seurat_clusters %in% c(1,4,6,7)), return.seurat = T)
seu_T_avg_nk <- AggregateExpression(seu_T %>% subset(seurat_clusters %in% c(8,10,11,12,13)), return.seurat = T)
seu_T_avg_cd4 <- AggregateExpression(seu_T %>% subset(seurat_clusters %in% c(0,2,3)), return.seurat = T)

# Visualize it as a heatmap with the gene_list you provide, change the object to visualize different ones
DoHeatmap(seu_T_avg_nk, features = gene_list, draw.lines = FALSE)  +  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 5, name = "RdBu"))) + theme_pubr(base_size = 10) + guides(colour = "none")

# Also visualize signature module scores from others, starting with Kyle Romine's
kyle_markerGenes <- read.csv(file = "Signatures/KR_TCellTypeSignatures.csv", header = T)

# Add module scores 
for (n in names(kyle_markerGenes)) {
  print(n)
  #n <- "HSPC"
  seu_T <- AddModuleScore(object = seu_T, features = kyle_markerGenes[n], name = n)
  colnames(seu_T@meta.data) <- gsub(str_c(n, "1$"), str_c(n, "_Score"), colnames(seu_T@meta.data))
  }

# Visualize
FeaturePlot(seu_T, features = c("CD4_Naïve_Score", "CD56_dim_NK_Score", "CD8_Term_Eff_Score", "CD8_GZMK_Exh_Score", "CD8_EM_Score", "CD8_Naïve_Score", "NKT_Score", "CD4_CM_Score", "MAIT_Score", "Tregs_Score", "CD56_Bright_NK_Score"))
ggsave("Tsubset_Kyle_FeaturePlots.pdf", width = 10, height = 8) # this was later moved into FeaturePlots/

# Based on the plots generated above, add T cell cluster names
Tcell.cluster.ids <- c("CD4 Memory", "CD8 Memory", "CD4 Naïve", "Treg", "CD8 Effector", "CD4 Effector Memory", "γδ T", "CD8 Exhausted", "CD56 Dim NK", "CD8 Naïve", "NK T", "CD56 Dim NK", "CD56 Bright NK", "Adaptive NK")
names(Tcell.cluster.ids) <- levels(seu_T)
seu_T <- RenameIdents(seu_T, Tcell.cluster.ids)
seu_T$celltype <- Idents(seu_T)

# See the levels
levels(seu_T$celltype)
levels(seu_T)
head(seu_T@meta.data)

# See the new annotated UMAP
DimPlot(seu_T, reduction = "umap", repel = T, group.by = "celltype", label = T) + theme(aspect.ratio = 1)

# Merge cell type annotations
meta <- seu@meta.data
tnk_meta <- seu_T@meta.data

meta$celltype <- as.character(meta$celltype)
tnk_meta$celltype <- as.character(tnk_meta$celltype)
meta$cell <- rownames(meta)
tnk_meta$cell <- rownames(tnk_meta)
meta$celltype[meta$cell %in% tnk_meta$cell] <- tnk_meta$celltype
seu@meta.data <- meta

# Turn celltype into a factor and add levels 
seu@meta.data$celltype <- factor(seu@meta.data$celltype,levels = c("Progenitors", "Early Erythroids", "Mid Erythroids", "Late Erythroids", "Pro Monocytes", "Monocytes", "Non Classical Monocytes", "cDC",  "pDC", "Pro B cells", "Pre-B", "B cells", "Plasma cells", "CD4 Naïve", "CD4 Effector Memory", "CD4 Memory", "Treg", "CD8 Naïve", "CD8 Effector", "CD8 Memory", "CD8 Exhausted",  "γδ T", "NK T", "Adaptive NK", "CD56 Bright NK", "CD56 Dim NK", "Cycling T-NK cells", "UD1", "UD2", "UD3"))
seu@active.ident <- seu$celltype 
DimPlot(seu, reduction = "umap", repel = T, label = T) + theme(aspect.ratio = 1)

# Add T-NK cells umap coordinates to the merged object. First, extract UMAP coordinates for TNK cells
UMAP_TNK_1 <- seu_T@reductions$umap@cell.embeddings[, 1]
UMAP_TNK_2 <- seu_T@reductions$umap@cell.embeddings[, 2]

# Add UMAP coordinates back to the metadata of the original Seurat object
seu$UMAP_TNK_1 <- NA
seu$UMAP_TNK_2 <- NA

# Populate only the TNK cells with their UMAP coordinates
seu$UMAP_TNK_1[Cells(seu_T)] <- UMAP_TNK_1
seu$UMAP_TNK_2[Cells(seu_T)] <- UMAP_TNK_2
seu@meta.data %>% head

# Save space
seu@assays$RNA@layers$scale.data <- NULL

# Save the merged seu as a new object
saveRDS(seu, file = "~/250411_SubsettedAnnotatedSeuratObject.rds")



# COMPARE/TMP
seu_old <- readRDS("~/250127_seu_annotated_merged_no-scale.rds")
seu

# Basic features - OK
seu_old
seu

# New one has additional "survival" column
seu_old@meta.data %>% head
seu@meta.data %>% head
setdiff(colnames(seu_old@meta.data), colnames(seu@meta.data))
setdiff(colnames(seu@meta.data), colnames(seu_old@meta.data))

# Same cells
identical(colnames(seu_old), colnames(seu))

# Not same celltype?
identical(seu_old$celltype, seu$celltype)
tmp_tib <- tibble(old_cell = colnames(seu_old),
                  new_cell = colnames(seu),
                  old_celltype = seu_old$celltype,
                  new_celltype = seu$celltype)
tmp_tib %>% filter(old_cell != new_cell)
tmp_tib %>% filter(old_celltype != new_celltype)


all( levels(seu_old$celltype) == levels(seu$celltype) )



