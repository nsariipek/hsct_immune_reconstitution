# Nurefsan Sariipek

# Since our project revolves around the T cells we wanted subset them and annotate them more closely.

#Load the saved seurat object
seu_diet <- readRDS("~/seu_diet.rds")

#Subset only T cells from all metadata
Tcell_subset <- subset(x = seu_diet, subset = seurat_clusters %in% c(0,1,2,5,6,11))

#Don't forget to run FindVariables again for more accurate analysis
Tcell_subset <- FindVariableFeatures(Tcell_subset, selection.method = "vst", nfeatures = 2000, verbose = FALSE)

#Decide dimensions according to elbow plot below
ElbowPlot(Tcell_subset)

Tcell_subset <- FindNeighbors(Tcell_subset, dims = 1:20)
Tcell_subset <- FindClusters(Tcell_subset, resolution = 1.1)

#Run UMAP
Tcell_subset <- RunUMAP(Tcell_subset2, dims = 1:20) 
DimPlot(Tcell_subset, reduction = "umap", label = TRUE) + theme(aspect.ratio = 1)

#Save for the next time
saveRDS(Tcell_subset, file= "~/Tcellsubset.RDS")

#You can save as diet seurat object
tc_diet <- DietSeurat(Tcell_subset, dimreducs = names(Tcell_subset@reductions))
saveRDS(tc_diet, file = "~/Tcellsubset_diet.rds")

#T cell Features to distinguish populations
FeaturePlot(tc_diet, features = c("CD8A", "CD8B", "CD4", "NCAM1","IL10","TGFB","GATA3","TCF7","SELL","CCR7","SELL","TMIGD2","LEF1","CD28","CD27"))

FeaturePlot(tc_diet, features = c("CD8A", "CD8B", "CD4", "NCAM1","IFNG", "CCL3", "PDCD1", "TIGIT", "LAG3", "HAVCR2", "CTLA4"))

#Run Find Markers
Tcell_subset_markers <- FindAllMarkers(Tcell_subset2, min.pct = .3, logfc.threshold = .3)

#Convert to a tibble and export it as a csv
Tcell_subset_markers_tibble <- as_tibble(Tcell_subset_markers)
write.csv(Tcell_subset_markers_tibble, file = "~/Tcell_subset_markers.csv")

# Heatmap the genes in your interest, providing them as gene_list, and change this as many times as you want.
# For the genes down below, we have used various published datasets, which will be referenced in our paper.

# Melanoma paper
gene_list = c("CD3E","CD4","CD8A","SELL","CCR7","IL7R","CD28","FAS","CD27","ITGAE","ITGAL","ITGAM","ITGAX","PDCD1","TIGIT","HAVCR2","LAG3","CTLA4","VTCN1","CD244","KLRG1","TNFRSF14","BTLA","CD160","CD38","ENTPD1","NT5E","CD69", "IL2RA","ICOS","TNFRSF4","TNFRSF9","HLA-DRA","CD40LG","GZMA","GZMB","GZMH","GZMK","GZMM","PRF1","NKG7","GNLY","IFNG","FASLG","TNF","IL17A","IL2","LEF1","TCF7","EOMES","TBX21","PRDM1","TOX","GATA3","ID2","ID3","NR4A1","ZNF683","FOXP3","MKI67","TOP2A","TRGV9","TRDV2","KLRB1","KLRC3")

gene_list = c("CD3E","CD4","CD8A","SELL","CCR7","IL7R","CD28","FAS","CD27","ITGAE","ITGAL","ITGAM","ITGAX","PDCD1","TIGIT","HAVCR2","LAG3","CTLA4","VTCN1","CD244","KLRG1","TNFRSF14","BTLA","CD160","CD38","ENTPD1","NT5E","CD69", "IL2RA","ICOS","TNFRSF4","TNFRSF9","HLA-DRA","CD40LG","GZMA","GZMB","GZMH","GZMK","GZMM","PRF1","NKG7","GNLY")


#CD8
gene_list = c("CD8", "CD8A", "CD8B", "CD4", "NKG7", "GNLY," "CST7", "PRF1", "GZMK," "GZMH", "GZMA," "GZMB," "IFNG," "CCL3", "PDCD1", "TIGIT", "LAG3", "HAVCR2", "CTLA4", "TCF7", "LEF1", "SELL," "CD27", "CD28", "CD57", "S1PR1", "VIM," "GPR183", "CCR7", "IL7R", "CCL3", "CCL3L1", "CCL4L2", "CCL4")

#CD4 
gene_list = c("CD4", "NKG7", "GNLY," "CCL4", "CST7", "PRF1", "GZMK," "GZMH," "GZMA," "GZMB," "IFNG," "CCL3", "PDCD1", "TIGIT," "LAG3", "HAVCR2", "CTLA4", "FOXP3", "IKZF2", "IL2RA", "CCR10", "TNFRSF4", "TIMP1", "USP10", "CCR7", "TCF7", "LEF1", "SELL," "IL7R", "ANXA1", "CD40LG", "RORA")

#CD4
gene_list = c("CCR7","SELL","TMIGD2","LEF1","ANXA1","LGALS1","TIMP1","S100A11","ANXA2","KLBR1","CCL5","GZMA","GZMK")

#CD56
gene_list = c("NCAM", "CD8","CD8A","CD8B", "GZMK", "XCL1", "IL7R", "TCF7", "GPR183", "GZMB", "PRF1", "CX3CR1", "CD7", "FCER1G", "KLRB1","KLRC2", "CD3E", "PATL2", "ZBTB38")

#th1 and 2 
gene_list = c("KLRD1","IFNGR1","CXCR3","CXCR6","CCR1","CCR5","STAT1","STAT4","TBX21","TNF","LTA","IFNG","IL2","IL12RB1","IL18R1","TNFSF11","HAVCR","CXCR4","BATF","IL4","CCR4","GATA3","IL5","IL13","CCR8","IRF4","AREG","STAT6","HAVCR1","PTGDR2","IL17RB","IL33","IL1R1", "AHR","CSF2","KLRB1","BATF","IL17A", "CCR4", "MAF", "IL17AF","CCR6","NFKBIZ","IL17F","IL21R","IRF4","IL21","IL22","RORA","RORC","STAT3","TBX21","PRF1","GZMB","GZMA")

#AverageExpression will normalize and scale the data when you set the return.seurat = TRUE

tc_diet_avg <- AverageExpression(tc_diet, return.seurat = TRUE)

tc_diet_avg_cd8 <- AverageExpression(tc_diet %>% subset(seurat_clusters %in% c(1,6)), return.seurat = TRUE)

tc_diet_avg_nk <- AverageExpression(tc_diet %>% subset(seurat_clusters %in% c(8,11,13)), return.seurat = TRUE)

tc_diet_avg_cd4 <- AverageExpression(tc_diet %>% subset(seurat_clusters %in% c(0,3,4,7,10,14)), return.seurat = TRUE)

#Visualize it either as a heatmap or dot plot with the gene_list you provide

#Heatmap
DoHeatmap(tc_diet_avg, features = gene_list, draw.lines = FALSE)  +  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 5, name = "RdBu"))) + theme_pubr(base_size = 10) + guides(colour = "none")

#DotPlot
DotPlot(tc_diet %>% subset(seurat_clusters %in% c(1,2,5,6,9)), features = gene_list, dot.scale = 10) + theme(axis.text.x = element_text(angle=90,vjust = 0.5), axis.text = element_text(size=11), text = element_text(size=11)) + scale_color_distiller(palette = "RdBu")

#Rename T cell Identities
Tcell.cluster.ids <- c("CD4 Memory","CD8 Effector","CD8 Memory","CD4 Naïve","Treg","CD8 Naïve","CD8 Effector","CD4 Naïve","CD56 Dim NK cells","CD8 Terminally Exhausted","CD4 Memory","NK T cells","γδ T lymphocytes","CD56 Bright NK cells", "CD4 Naïve")
names(Tcell.cluster.ids) <- levels(tc_diet)
tc_diet <- RenameIdents(tc_diet, Tcell.cluster.ids)
tc_diet@meta.data$celltype = Idents(tc_diet)

# See the new annotated UMAP
table(seu_diet_merged)
Idents(tc_diet) = "seurat_clusters"
DimPlot(tc_diet, reduction = "umap") + theme(aspect.ratio = 1)


# Merge cell type annotation
meta = seu_diet@meta.data
tnk_meta = tc_diet@meta.data

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




