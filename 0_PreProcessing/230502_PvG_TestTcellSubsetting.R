# Peter van Galen, 230503
# Check T cell subsetting and dimensionality reduction from Nurefsan's data


library(tidyverse)
library(Seurat)
library(cowplot)

# Make the objects more manageable ----------------------------------------------------------------

# Run this section only once, use diet objects next time

# First, copy data to persisent disk (in RStudio Terminal tab)
#gsutil cp gs://fc-3783b423-62ac-4c69-8c2f-98cb0ee4503b/p53_Nurefsan.rds .
#gsutil cp gs://fc-3783b423-62ac-4c69-8c2f-98cb0ee4503b/Tcellsubset.RDS .

seu <- readRDS("~/p53_Nurefsan.rds")
seu_diet <- DietSeurat(seu, dimreducs = names(seu@reductions))
saveRDS(seu_diet, file = "~/p53_Nurefsan_diet.rds")

tc <- readRDS("~/Tcellsubset.RDS")
tc_diet <- DietSeurat(tc, dimreducs = names(tc@reductions))
saveRDS(tc_diet, file = "~/Tcellsubset_diet.rds")


# Check variable genes ----------------------------------------------------------------------------

seu <- readRDS("~/p53_Nurefsan_diet.rds")
tc <- readRDS("~/Tcellsubset_diet.rds")

LabelPoints(plot = VariableFeaturePlot(seu), points = head(VariableFeatures(seu), 20), repel = T, xnudge = 0, ynudge = 0) +
  theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5)) + ggtitle("Variable genes")
LabelPoints(plot = VariableFeaturePlot(tc), points = head(VariableFeatures(tc), 20), repel = T, xnudge = 0, ynudge = 0) +
  theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5)) + ggtitle("Variable genes")

# Some plots
DimPlot(tc) + theme(aspect.ratio = 1)
p1 <- FeaturePlot(tc, features = "NCAM1") + theme(aspect.ratio = 1) + theme(aspect.ratio = 1)
p2 <- FeaturePlot(tc, features = "CD8A") + theme(aspect.ratio = 1) + theme(aspect.ratio = 1)
p3 <- FeaturePlot(tc, features = "CD4") + theme(aspect.ratio = 1) + theme(aspect.ratio = 1)
plot_grid(p1,p2,p3)

# Check variable features
VariableFeatures(tc)
tc2 <- tc
tc2 <- FindVariableFeatures(tc2, nfeatures = 2000, verbose = FALSE)
all( VariableFeatures(tc2) == VariableFeatures(tc) )






# EVERYTHING BELOW IS NUREFSANS UNALTERED CODE (p53.R) --------------------------------------------

#Nurefsan Sariipek, 230405
#Analysis of data of TP53 project
#Load libraries
library(tidyverse)
library(Seurat)
library(ggplot2)
#library(harmony)
#library(randomcoloR)
library(readxl)
library(data.table)
#library(ggforce)
#library(cloudml)
# Start with a clean slate
rm(list=ls())
# Parameters to interact with Google bucket
project <- Sys.getenv('WORKSPACE_NAMESPACE')
workspace <- Sys.getenv('WORKSPACE_NAME')
bucket <- Sys.getenv('WORKSPACE_BUCKET')
data_dir <- gs_data_dir_local(bucket)

#Load matrices of all samples
Samples <- c("1013_MNC", "1195_MNC", "1285_MNC", "1347_CD3", "1347_MNC", "1677_CD3", "1677_MNC", "1732_MNC", "1811_MNC", "1811_CD3", "1953_MNC", "1953_CD3", "1972_CD3", "1972_MNC", "2220_MNC", "2379_CD3", "2379_MNC", "2434_MNC", "2446_MNC", "2518_CD3", "2518_MNC", "2599_CD3", "2599_MNC", "2621_CD3", "2621_MNC", "2645_MNC", "2737_CD3", "2737_MNC", "4618_MNC","5641_MNC", "6174_MNC", "6174_CD3", "6244_MNC", "6244_CD3", "9185_MNC", "9185_CD3", "9596_CD3", "9596_MNC", "9355_MNC", "9931_MNC", "9931_CD3", "25802_MNC", "25802_CD3", "25809_MNC")

matrices_ls <- lapply(Samples, function(x) Read10X((data.dir = paste0("gs/fc-3783b423-62ac-4c69-8c2f-98cb0ee4503b/", x, "/sample_filtered_feature_bc_matrix/"))))

#Merge Files and Create Seurat Object
seu_ls <- lapply(matrices_ls, function(x) CreateSeuratObject(x))
for (n in 1:length(Samples)){
  object <- seu_ls[[n]]
  object@meta.data$orig.ident <- Samples[n]
  seu_ls[[n]] <- object }
seu <- merge(seu_ls[[1]], seu_ls[2:length(seu_ls)])
#use this code the clean the memory
gc()

# Add QC Metrics
seu <- PercentageFeatureSet(seu, "^MT-", col.name = "percent_mito")
seu <- PercentageFeatureSet(seu, "^RB[SL]", col.name = "percent_ribo")
#seu <- PercentageFeatureSet(seu, "^HB[^(P)]", col.name = "percent_hb") 

#Count the cells before filtering 
table(seu$orig.ident)

#QC Control Thresholds
seu <- subset(seu, subset = nFeature_RNA > 250 & nCount_RNA > 500 & percent_mito < 20)   

#Count the cells after filtering
table(seu$orig.ident)

#Add a column with cohort identities
seu $cohort <- case_when(grepl("9596|2737|2379|2434|2518|4618|6174|9931|1953|9355|1013|25809", seu$orig.ident) ~ "cohort2",
                         grepl("2446|25802|2645|1972|2220|2621|9185|2599", seu$orig.ident) ~ "cohort1",
                         grepl("1677|5641|1732|1811|1195|1347|1285|6244", seu$orig.ident) ~ "cohort3" )


#Add timepoints of samples as a column to metadata column 
seu $status <- case_when(grepl("9596|2379|4618|9355|2446|1972|1677|1195|5641", seu$orig.ident) ~ "pre_transplant",
                         grepl("25809|2434|2518|6174|9931|1013|25802|2645|2220|2621|9185|2599|1732|1285|6244", seu$orig.ident) ~ "remission",
                         grepl("2737|1953|1811|1347", seu$orig.ident) ~ "relapse")



#use this line to convert to factor
seu$cohort <- as.factor(seut@meta.data$cohort)
seu$status <- as.factor(seu@meta.data$status)

#use this library to visualize the numbers in different cohorts
#library("janitor")
as_tibble(seu@meta.data) %>% tabyl(status)
as_tibble(seu@meta.data) %>% tabyl(cohort)
table(seu$cohort)
table(seu$cohort)
seu@meta.data

# Visualize QC metrics as a violin plot
feats <- c("nFeature_RNA", "nCount_RNA", "percent_mito", "percent_ribo")
VlnPlot(seu, group.by = "orig.ident", features = feats, pt.size = 0.1, ncol = 3) + NoLegend()
FeatureScatter(seu, "nCount_RNA", "nFeature_RNA", group.by = "orig.ident", pt.size = 0.5)
gc()

#Data Normalization
seu <- NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 10000)
gc()

# Identify the most variable genes
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
top10 <- head(VariableFeatures(seu), 10)
write.csv(top10, file = "top10.csv")
plot1 <- VariableFeaturePlot(seu)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 
all.genes <- rownames(seu)
gc()
seu <- ScaleData(seu, features = all.genes)
gc()
seu <- RunPCA(seu, features = VariableFeatures(object = seu))
gc()

#Determine Dimensionality of Data
ElbowPlot(seu)
VizDimLoadings(seu, dims = 1:2, reduction = "pca")

#Visualize PCA results a few different ways
DimPlot(seu, reduction = "pca")
DimHeatmap(seu, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(seu, dims = 1:15, cells = 500, balanced = TRUE)
Idents(seu) = "seurat_clusters"

#Determine the dimensionality of the dataset
seu <- JackStraw(seu, num.replicate = 100)
seu <- ScoreJackStraw(seu, dims = 1:20)
JackStrawPlot(seu, dims = 1:20)

seu <- FindNeighbors(seu, dims = 1:20)
gc()

seu <- FindClusters(seu, resolution = 1)
ahead(Idents(seu), 5)
gc() 

#Run UMAP
seu <- RunUMAP(seu, dims = 1:25)

#Visualize UMAP
DimPlot(seu, reduction = "umap", label = TRUE) + theme(aspect.ratio = 1)

#Visualize UMAPs with different identities
UMAP_sample <- DimPlot(seu, reduction = "umap", group.by = "orig.ident") + theme(aspect.ratio = 1)
UMAP_sample

UMAP_cohort <- DimPlot(seu, reduction = "umap", group.by = "cohort") + theme(aspect.ratio = 1)
UMAP_cohort

UMAP_status <- DimPlot(seu, reduction = "umap", group.by = "status") + theme(aspect.ratio = 1)
UMAP_status
#split the UMAPs
UMAP2 <- DimPlot(seu, reduction = "umap", group.by = "orig.ident", split.by = "cohort") + theme(aspect.ratio = 1)
UMAP2
UMAP3 <- DimPlot(seu, reduction = "umap", group.by = "orig.ident", split.by = "status") + theme(aspect.ratio = 1)
UMAP3

gc()

#saveRDS as an object at this point and do not modify it after this part meaning do not re-save it
saveRDS(seu, file = "~/p53_Nurefsan.rds")

#next time you can just load your RDS assigning to a object for needed analysis
seu <-readRDS("~/p53_Nurefsan.rds")

gc()

#Run FindAllMarkers for to see msot expressed genes
seu_markers <- FindAllMarkers(seu, min.pct = .3, logfc.threshold = .3)

#Save as tibble and export as a csv for next time and also to export and work on cell type annotations
seu_markers_tib <- as_tibble(seu_markers)
write.csv(seu_markers_tib, file = "~/seu_markers_tib.csv")

#Load Marker Genes that you have saved before
seu_markers <- read.csv(file = "~/seu_markers_tib.csv")

#To visualize and help with annotating clusters run the lines below
#General Features
FeaturePlot(seu, features = c("CD34","MPO", "CD14", "MS4A1", "CD3G","CD8B"))

#Feature Plots
#T cell Features
FeaturePlot(seu, features = c("CD8A", "TCF7", "TOX", "HAVCR2", "CXCR3",
                              "SLAMF6", "CD3E", "CD4", "SELL", "CD44", "PDCD1",
                              "FOXP3", "GZMB", "GZMK", "LAG3", "CD101",
                              "CXCR5", "KLRG1", "IFNG", "TNF")) 
#Myeloid and B cell Features
FeaturePlot(seu, features = c("CCR2", "CD14", "CD33", "CD34", "THY1", "IL3RA",
                              "CX3CR1", "MME", "PTPRC", "ITGAX", "CD80", "CD19",
                              "ITGAL", "XCR1", "LY6C", "CSF1R", "ADGRE1")) 

#HSC markers
FeaturePlot(seu, features = c("VIM", "FLT3", "CD34", "ITGAL", "THY1",
                              "PTPRC", "KIT", "SLAMF1", "MME", "SLAMF2", "MPO")) 
#B cell markers
FeaturePlot(seu, features = c("CD19", "MS4A1", "PDCD1LG2", "NT5E", "FCER2", "SDC1",
                              "PAX5", "TCF3", "CD80", "SPIB", "BCL6")) 

#Compare Groups 
seu@active.ident <- seu$seurat_clusters  #Set ident to desired comparison group
#2 clusters potentially Monocytes
C20vC21 <- FindMarkers(seu, ident.1 = "20", 
                     ident.2 = "21",
                     logfc.threshold = 0.25,
                     test.use = "wilcox",
                     min.pct = 0.1)
#2 clusters that might be pro-B cells
C24vC12 <- FindMarkers(seu, ident.1 = "20", 
                       ident.2 = "21",
                       logfc.threshold = 0.25,
                       test.use = "wilcox",
                       min.pct = 0.1)

#Rename Clusters
seu.cluster.ids <- c("T cells","Monocytes","Blasts","T cells","T cells","Late Erythroids","Monocytes","Mid Erythroids", "B cells","T cells","Pro Monocytes", "Early Erythroids","Early Erythroids","Late Erythroids","Blasts","GMP","Non Classical Monocytes", "B cells", "Monocytes", "Monocytes", "cDC","MEP","B cells", "HSPCs", "Pro B cells", "Plasma Cells", "Blasts","Blasts","B cells", "Late Erythroids", "pDC", "Late Erythroids", "Doublets")
names(seu.cluster.ids) <- levels(seu)
seu <- RenameIdents(seu, seu.cluster.ids)
seu@meta.data$celltype = Idents(seu)

#Make new heatmap with idents
mycolors <- distinctColorPalette(k = 34)
pie(rep(1, 34), col = mycolors) 
DimPlot(seu, reduction = "umap", repel = T, group.by = "celltype", cols = mycolors, label = T) + theme(aspect.ratio = 1)
UMAP_cohort2 <- DimPlot(seu, reduction = "umap", group.by = "cohort") + theme(aspect.ratio = 1)
UMAP_cohort2
UMAP_status2 <- DimPlot(seu, reduction = "umap", group.by = "status") + theme(aspect.ratio = 1)
UMAP_status2

########## Subset T cells  ########
Tcell_subset2 <- subset(x = seu, subset = seurat_clusters %in% c(0,1,2,5,6,11))
Tcell_subset2 <- FindVariableFeatures(Tcell_subset2, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
ElbowPlot(Tcell_subset2)
Tcell_subset2 <- FindNeighbors(Tcell_subset2, dims = 1:20)
Tcell_subset2 <- FindClusters(Tcell_subset2, resolution = 1.1)
gc()

#Run UMAP
Tcell_subset2 <- RunUMAP(Tcell_subset2, dims = 1:20) 

DimPlot(Tcell_subset2, reduction = "umap", label = TRUE) + theme(aspect.ratio = 1)

saveRDS(Tcell_subset2, file= "~/Tcellsubset.RDS")

Tcell_subset <- readRDS("~/Tcellsubset.RDS")

# T cell Features
FeaturePlot(Tcell_subset, features = c("CD8A", "CD8B", "TCF7", "TOX", "HAVCR2", "CXCR3",
                                       "SLAMF6", "CD3E", "CD4", "SELL", "CD44", "PDCD1",
                                       "FOXP3", "GZMB", "GZMK", "LAG3", "CD101",
                                       "CXCR5", "KLRG1", "IFNG", "TNF")) 
#free the memory before running this line
gc()

#Run Find Markers
Tcell_subset_markers <- FindAllMarkers(Tcell_subset2, min.pct = .3, logfc.threshold = .3)
#Convert to a tibble and export it as a csv
Tcell_subset_markers_tibble <- as_tibble(Tcell_subset_markers)
write.csv(Tcell_subset_markers_tibble, file = "~/Tcell_subset_markers.csv")


Tcell_UMAP_sample <- DimPlot(Tcell_subset, reduction = "umap", group.by = "orig.ident") + theme(aspect.ratio = 1)
Tcell_UMAP_sample
Tcell_UMAP_status <- DimPlot(Tcell_subset, reduction = "umap", group.by = "status") + theme(aspect.ratio = 1)
Tcell_UMAP_status
Tcell_UMAP_cohort <- DimPlot(Tcell_subset, reduction = "umap", group.by = "cohort") + theme(aspect.ratio = 1) 
Tcell_UMAP_cohort


#For annotation purposes use Kyle's signatures 
Kylesgenes <- read.csv(file = "~/KR_CellTypeSignatures.csv", header = T)
KyleTcells <- read.csv(file = "~/KR_TCellTypeSignatures.csv", header = T)
#Add Module Scores for cell type annotation
for (n in names(KyleTcells)) {
  print(n)
  #n <- "HSPC"
Tcell_subset <- AddModuleScore(object = Tcell_subset, features = KyleTcells[n], name = n)
  colnames(Tcell_subset@meta.data) <- gsub(str_c(n, "1$"), str_c(n, "_Score"), colnames(Tcell_subset@meta.data))
}

FeaturePlot(Tcell_subset, features = c("CD4_Naïve_Score", "CD56_dim_NK_Score", "CD8_Term_Eff_Score", "CD8_GZMK_Exh_Score", "CD8_EM_Score", "CD8_Naïve_Score", "NKT_Score", "CD4_CM_Score", "MAIT_Score", "Tregs_Score", "CD56_Bright_NK_Score"))



#CD4_CM
FeaturePlot(seu, features = c("LTB","AL136456.1","IL7R", "FRY","AC079793.1","ADAM19", "PTPN13","BCL2","VIM","RORA-AS1","AQP3","LINC00513","GAS5", "INPP4B", "LTB","ANK3","SERINC5", "ADAM19", "SLC2A3", "IL7R", "CDC14A","FAAH2", "ARHGAP15", "PAG1", "SAMSN1","ICOS"))
#cd4_naive  
FeaturePlot(Tcell_subset, features = c("INPP4B","LTB","ANK3","SERINC5","ADAM19","SLC2A3","IL7R","CDC14A","FAAH2","ARHGAP15","PAG1","SAMSN1",
"ICOS"))
#cd8 naive t cells
FeaturePlot(Tcell_subset, features = c("NELL2","THEMIS","SIK3","LINC01619","GZMK","ATXN1","IL21R","PARP8","WWOX","CBLB","GGA2","PDE3B","P2RY8",
"PPP1R16B"))

#Rename Tcell Identities
Tcell.cluster.ids <- c("",  )
names(Tcell.cluster.ids) <- levels(Tcell_subset)
Tcell_subset <- RenameIdents(Tcell_subset, Tcell.cluster.ids)
Tcell_subset@meta.data$celltype = Idents(Tcell_subset)

#subset the blasts populations
blasts <- subset(seu, subset =celltype %in% c("Blasts"))

#Compare differential gene expression between groups
seu@active.ident <- seu$status
Idents(blasts) <- blasts$status
Analysis1<- FindMarkers(blasts, ident.1 = "relapse" , 
                     ident.2 = "remission",
                     logfc.threshold = 0.25,
                     test.use = "wilcox",
                     min.pct = 0.1)
Analysis1_tibble <- as_tibble(Analysis1)
Analysis1_tibble = Analysis1 %>% rownames_to_column("gene")
write.csv(Analysis1_tibble, file = "~/Analysis1.csv", quote = F, append = F, row.names = F)


seu@active.ident <- seu$cohort
Idents(blasts) <- blasts$cohort
Analysis4<- FindMarkers(blasts, ident.1 = "remission" , 
                        ident.2 = "relapse",
                        logfc.threshold = 0.25,
                        test.use = "wilcox",
                        min.pct = 0.1)
Analysis4_tibble <- as_tibble(Analysis4)
Analysis4_tibble = Analysis4 %>% rownames_to_column("gene")
write.csv(Analysis4_tibble, file = "~/Analysis4.csv", quote = F, append = F, row.names = F)



seu@active.ident <- seu$cohort
Idents(Tcell_subset) <- Tcell_subset$cohort
Analysis7<- FindMarkers(Tcell_subset, ident.1 = "relapse" , 
                                ident.2 = "remission",
                                logfc.threshold = 0.25,
                                test.use = "wilcox",
                                min.pct = 0.1)
Analysis7_tibble <- as_tibble(Analysis7)
Analysis7_tibble = Analysis7 %>% rownames_to_column("gene")
write.csv(Analysis7_tibble, file = "~/Analysis7.csv", quote = F, append = F, row.names = F)


seu@active.ident <- Tcell_subset2$cohort
Idents(Tcell_subset2) <- Tcell_subset2$cohort
Analysis5<- FindMarkers(Tcell_subset2, ident.1 = "relapse" , 
                        ident.2 = "remission",
                        logfc.threshold = 0.25,
                        test.use = "wilcox",
                        min.pct = 0.1)
Analysis5_tibble <- as_tibble(Analysis5)
Analysis5_tibble = Analysis5 %>% rownames_to_column("gene")
write.csv(Analysis5_tibble, file = "~/Analysis5.csv", quote = F, append = F, row.names = F)

seu@active.ident <- Tcell_subset2$status
Idents(Tcell_subset2) <- Tcell_subset2$status
Analysis6<- FindMarkers(Tcell_subset2, ident.1 = "relapse" , 
                        ident.2 = "remission",
                        logfc.threshold = 0.25,
                        test.use = "wilcox",
                        min.pct = 0.1)
Analysis6_tibble <- as_tibble(Analysis6)
Analysis6_tibble = Analysis6 %>% rownames_to_column("gene")
write.csv(Analysis6_tibble, file = "~/Analysis6.csv", quote = F, append = F, row.names = F)

