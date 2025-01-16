# Nurefsan Sariipek
# Annotating T cells in a more granular way
library(tidyverse)
library(Seurat)
library(janitor)

# Start with a clean slate
rm(list=ls())

# Set working directory 
setwd("~/TP53_ImmuneEscape/2_Annotate/")

# Load the saved Seurat object from the last step that contains only 20% of the cells
seu <- readRDS("~/250113_SplittedSeuratObject.rds")

# Subset only T cells from all metadata
Tcell_subset <- subset(x = seu, subset = seurat_clusters %in% c(0,1,3,5,11,24,30))

# Run FindVariables again (replace the ones inherited from the complete object)
Tcell_subset <- FindVariableFeatures(Tcell_subset, selection.method = "vst", nfeatures = 2000, verbose = FALSE)

# Decide dimensions according to elbow plot below
ElbowPlot(Tcell_subset)
Tcell_subset <- FindNeighbors(Tcell_subset, dims = 1:20)
Tcell_subset <- FindClusters(Tcell_subset, resolution = 1.1)

# Run UMAP
Tcell_subset <- RunUMAP(Tcell_subset, dims = 1:20)
#Visualize the UMAP
umap_plot <- DimPlot(Tcell_subset, reduction = "umap", label = TRUE) + theme(aspect.ratio = 1)
pdf("UMAP_Tcells.pdf", width = 10, height = 8)  # Adjust width and height as needed
print(umap_plot)
dev.off()

# Save for the next time
saveRDS(Tcell_subset, file= "~/Tcellsubset.RDS")

# Load the saved Seurat object
Tcell_subset <- readRDS("~/250116_Tcellsubset.RDS")

# Visualize the numbers in different cohorts
as_tibble(Tcell_subset@meta.data) %>% tabyl(sample_status)
as_tibble(Tcell_subset@meta.data) %>% tabyl(cohort)
as_tibble(Tcell_subset@meta.data) %>% tabyl(library_type)

# Run Find Markers
Tcell_subset_markers <- FindAllMarkers(Tcell_subset, min.pct = .3, logfc.threshold = .3)

# Convert to a tibble and export it as a csv
Tcell_subset_markers_tibble <- as_tibble(Tcell_subset_markers)
write.csv(Tcell_subset_markers_tibble, file = "~/Tcell_subset_markers.csv")


# T cell Features to distinguish populations
FeaturePlot(Tcell_subset, features = c("CD8A", "CD8B", "CD4", "NCAM1","IL10","TGFB","GATA3","TCF7","SELL","CCR7","SELL","TMIGD2","LEF1","CD28","CD27"))

FeaturePlot(Tcell_subset, features = c("CD8A", "CD8B", "CD4", "NCAM1","IFNG", "CCL3", "PDCD1", "TIGIT", "LAG3", "HAVCR2", "CTLA4"))

#Naive
FeaturePlot(Tcell_subset, features = c("CCR7","TCF7","LEF1", "SELL"))

FeaturePlot(Tcell_subset, features = c("CD3E"))

#CD8 Naive 
FeaturePlot(Tcell_subset, features = c("CD8A", "CD8B", "LEF1", "SELL", "CCR7", "TCF7", "CD27", "CD28","S1PR1"))

#Effector 
FeaturePlot(Tcell_subset, features = c("GZMA","GZMB","PRF1","CX3CR1"))

#Cytotoxicity 
FeaturePlot(Tcell_subset, features = c("CX3CR1","GMZH","GNLY","FGFBP2"))

#Cytotoxic
FeaturePlot(Tcell_subset, features = c("FCGR3A", "KLRG1", "PRF1", "GZMB"))

#TIL Atlas
#cd8_tpex
FeaturePlot(Tcell_subset, features =c("LAG3", "XCL1", "CRTAM","IFNG", "CCL4", "PDCD1", "DUSP4", "CD8A", "ZEB2", "NR4A2", "SLA", "NKG7", "TIGIT", "CTSW", "TNFRSF9", "TOX", "LYST", "TNFSF4", "CCL3", "GZMB", "RAB27A", "PRF1", "CD70", "PLSCR1","CXCL13"))

#CD4 naive like
FeaturePlot(Tcell_subset, features = c("CCR7", "SELL", "CD40LG","IL7R", "TCF7", "LEF1", "GPR183", "KLRB1", "LTB", "MAL", "PASK", "AQP3", "TRAT1"))      

#Cytotoxic
FeaturePlot(Tcell_subset, features = c("NKG7", "CCL4", "CST7", "GZMA", "GZMB", "IFNG", "CCL3"))

#Effector memory
FeaturePlot(Tcell_subset, features = c("KLRG1", "LYAR", "GZMK", "GZMM", "TXNIP","FCRL6"))

#CD8, Memory cells 
FeaturePlot(Tcell_subset, features = c("TCF7","CCR7","IL7R","CD8A","CD8B","SELL","GZMA"))

#Memory 
FeaturePlot(Tcell_subset, features = c("TCF7","CCR7","IL7R","SELL"))

#T-NK
FeaturePlot(Tcell_subset, features = c("TCF7","GZMM","KLRC3","KLRB1","GNLY"))

#Exhausted
FeaturePlot(Tcell_subset, features = c("LAYN", "ITGAE", "PDCD1", "CTLA4", "HAVCR2", "LAG3", "TIGIT"))

#Exhaustion
FeaturePlot(Tcell_subset, features = c("PDCD1", "TIGIT", "HAVCR2", "LAG3"))

#Terminally exhausted
FeaturePlot(Tcell_subset, features = c("TRAV24","TRBV5-5","CCDC50","HLA-DRB1", "VCAM1", "CTSW", "NKG7", "CXCL13", "HLA-DPA1", "DUSP4", "GZMB","HLA-DRA"))

#CD8 T Effector memory cells
FeaturePlot(Tcell_subset, features = c("S1PR1", "GPR183", "CCR7", "ANXA1", "TCF7", "IL7R", "MBP", "VIM","GYG1","GPR183"))

#CD8 T progenitors
FeaturePlot(Tcell_subset, features = c("CAV1", "GNG4", "XCL1", "CRTAM", "CXCL13", "GEM","XCL2"))
#TREG
FeaturePlot(Tcell_subset, features = c("FOXP3","IL2RA","SELL","TNFRSF4"))
#MAIT cells
FeaturePlot(Tcell_subset, features = c("KLRB1", "NCR3", "ZBTB16", "SLC4A10", "RORC"))

# Heatmap the genes in your interest, providing them as gene_list, and change this as many times as you want.
# For the genes down below, we have used various published datasets, which will be referenced in our paper.

# Melanoma paper
gene_list = c("CD3E","CD4","CD8A","SELL","CCR7","IL7R","CD28","FAS","CD27","ITGAE","ITGAL","ITGAM","ITGAX","PDCD1","TIGIT","HAVCR2","LAG3","CTLA4","VTCN1","CD244","KLRG1","TNFRSF14","BTLA","CD160","CD38","ENTPD1","NT5E","CD69", "IL2RA","ICOS","TNFRSF4","TNFRSF9","HLA-DRA","CD40LG","GZMA","GZMB","GZMH","GZMK","GZMM","PRF1","NKG7","GNLY","IFNG","FASLG","TNF","IL17A","IL2","LEF1","TCF7","EOMES","TBX21","PRDM1","TOX","GATA3","ID2","ID3","NR4A1","ZNF683","FOXP3","MKI67","TOP2A","TRGV9","TRDV2","KLRB1","KLRC3")

gene_list = c("CD3E","CD4","CD8A","SELL","CCR7","IL7R","CD28","FAS","CD27","ITGAE","ITGAL","ITGAM","ITGAX","PDCD1","TIGIT","HAVCR2","LAG3","CTLA4","VTCN1","CD244","KLRG1","TNFRSF14","BTLA","CD160","CD38","ENTPD1","NT5E","CD69", "IL2RA","ICOS","TNFRSF4","TNFRSF9","HLA-DRA","CD40LG","GZMA","GZMB","GZMH","GZMK","GZMM","PRF1","NKG7","GNLY")


#CD8 genes
gene_list = c("CD8", "CD8A", "CD8B", "CD4", "NKG7", "GNLY", "CST7", "PRF1", "GZMK", "GZMH", "GZMA", "GZMB", "IFNG" ,"CCL3", "PDCD1", "TIGIT", "LAG3", "HAVCR2", "CTLA4", "TCF7", "LEF1", "SELL", "CD27", "CD28", "CD57", "S1PR1", "VIM", "GPR183", "CCR7", "IL7R", "CCL3", "CCL3L1", "CCL4L2", "CCL4")

#CD4 genes
gene_list = c("CCR7","SELL","TMIGD2","LEF1","ANXA1","LGALS1","TIMP1","S100A11","ANXA2","KLBR1","CCL5","GZMA","GZMK")

#CD56 NK cells genes
gene_list = c("NCAM", "CD8","CD8A","CD8B", "GZMK", "XCL1", "IL7R", "TCF7", "GPR183", "GZMB", "PRF1", "CX3CR1", "CD7", "FCER1G", "KLRB1","KLRC2", "CD3E", "PATL2", "ZBTB38")

#Th1 and Th2 
gene_list = c("KLRD1","IFNGR1","CXCR3","CXCR6","CCR1","CCR5","STAT1","STAT4","TBX21","TNF","LTA","IFNG","IL2","IL12RB1","IL18R1","TNFSF11","HAVCR","CXCR4","BATF","IL4","CCR4","GATA3","IL5","IL13","CCR8","IRF4","AREG","STAT6","HAVCR1","PTGDR2","IL17RB","IL33","IL1R1", "AHR","CSF2","KLRB1","BATF","IL17A", "CCR4", "MAF", "IL17AF","CCR6","NFKBIZ","IL17F","IL21R","IRF4","IL21","IL22","RORA","RORC","STAT3","TBX21","PRF1","GZMB","GZMA")

# AverageExpression will normalize and scale the data when you set the return.seurat = TRUE
Tcell_subset_avg <- AverageExpression(Tcell_subset, return.seurat = TRUE)

Tcell_subset_avg_cd8 <- AverageExpression(Tcell_subset %>% subset(seurat_clusters %in% c(1,6)), return.seurat = TRUE)

Tcell_subset_avg_nk <- AverageExpression(Tcell_subset %>% subset(seurat_clusters %in% c(8,11,13)), return.seurat = TRUE)

Tcell_subset_avg_cd4 <- AverageExpression(Tcell_subset %>% subset(seurat_clusters %in% c(0,3,4,7,10,14)), return.seurat = TRUE)

# Visualize it either as a heatmap or dot plot with the gene_list you provide, change the object to visualize different ones
# Heatmap
DoHeatmap(Tcell_subset_avg, features = gene_list, draw.lines = FALSE)  +  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 5, name = "RdBu"))) + theme_pubr(base_size = 10) + guides(colour = "none")

# DotPlot
DotPlot(Tcell_subset %>% subset(seurat_clusters %in% c(1,2,5,6,9)), features = gene_list, dot.scale = 10) + theme(axis.text.x = element_text(angle=90,vjust = 0.5), axis.text = element_text(size=11), text = element_text(size=11)) + scale_color_distiller(palette = "RdBu")

# Rename T cell Identities
Tcell.cluster.ids <- c("CD4 Memory","CD8 Effector","CD8 Memory","CD4 Naïve","Treg","CD8 Naïve","CD8 Effector","CD4 Naïve","CD56 Dim NK cells","CD8 Terminally Exhausted","CD4 Memory","NK T cells","γδ T lymphocytes","CD56 Bright NK cells", "CD4 Naïve")
names(Tcell.cluster.ids) <- levels(Tcell_subset)
Tcell_subset <- RenameIdents(Tcell_subset, Tcell.cluster.ids)
Tcell_subset@meta.data$celltype = Idents(Tcell_subset)

# See the new annotated UMAP
table(seu_diet_merged)
Idents(Tcell_subset) = "seurat_clusters"
DimPlot(Tcell_subset, reduction = "umap") + theme(aspect.ratio = 1)

# Merge cell type annotations
meta = seu_diet@meta.data
tnk_meta = Tcell_subset@meta.data

meta$celltype = as.character(meta$celltype)
tnk_meta$celltype = as.character(tnk_meta$celltype)
meta$cell = rownames(meta)
tnk_meta$cell = rownames(tnk_meta)
meta$celltype[meta$cell %in% tnk_meta$cell] = tnk_meta$celltype
seu_diet@meta.data = meta

# Save the merged seu as a new object
seu_diet_merged <- DietSeurat(seu_diet, dimreducs = names(seu_diet@reductions))
saveRDS(seu_diet_merged, file = "~/seu_diet_merged.rds")
seu_diet_merged <- readRDS("~/seu_diet_merged.rds")