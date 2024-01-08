# Nurefsan Sariipek, 230405
# Analysis of data of TP53 project
# Load the needed libraries
library(tidyverse)
library(Seurat)
library(ggplot2)
library(randomcoloR)
library(readxl)
library(data.table)
library(janitor)

# Set the home directory
setwd("~/")

# Load matrices of all samples
Samples <- c("1013_MNC", "1195_MNC", "1285_MNC", "1347_CD3", "1347_MNC", "1677_CD3", "1677_MNC", "1732_MNC", "1811_MNC", "1811_CD3", "1953_MNC", "1953_CD3", "1972_CD3", "1972_MNC", "2220_MNC", "2379_CD3", "2379_MNC", "2434_MNC", "2446_MNC", "2518_CD3", "2518_MNC", "2599_CD3", "2599_MNC", "2621_CD3", "2621_MNC", "2645_MNC", "2737_CD3", "2737_MNC", "4618_MNC","5641_MNC", "6174_MNC", "6174_CD3", "6244_MNC", "6244_CD3", "9185_MNC", "9185_CD3", "9596_CD3", "9596_MNC", "9355_MNC", "9931_MNC", "9931_CD3", "25802_MNC", "25802_CD3", "25809_MNC")

matrices_ls <- lapply(Samples, function(x) Read10X((data.dir = paste0("gs/fc-3783b423-62ac-4c69-8c2f-98cb0ee4503b/", x, "/sample_filtered_feature_bc_matrix/"))))

#Merge Files and Create Seurat Object
seu_ls <- lapply(matrices_ls, function(x) CreateSeuratObject(x))
for (n in 1:length(Samples)){
  object <- seu_ls[[n]]
  object@meta.data$orig.ident <- Samples[n]
  seu_ls[[n]] <- object }
seu <- merge(seu_ls[[1]], seu_ls[2:length(seu_ls)])   

# Add QC Metrics
seu <- PercentageFeatureSet(seu, "^MT-", col.name = "percent_mito")
seu <- PercentageFeatureSet(seu, "^RB[SL]", col.name = "percent_ribo")

# Count the cells before filtering 
table(seu$orig.ident)

# QC Control Thresholds
seu <- subset(seu, subset = nFeature_RNA > 250 & nCount_RNA > 500 & percent_mito < 20)   

# Count the cells after filtering
table(seu$orig.ident)

# Visualize QC metrics as a violin plot
                 
feats <- c("nFeature_RNA", "nCount_RNA", "percent_mito", "percent_ribo")
VlnPlot(seu_diet, group.by = "patient_identity", features = feats, pt.size = 0.1, ncol = 3) + NoLegend()

FeatureScatter(seu_diet_merged, "nCount_RNA", "nFeature_RNA", group.by = "id", pt.size = 0.5) + theme(aspect.ratio = 1) +
  geom_vline(xintercept=250, col="black") +
  geom_hline(yintercept=500, col="black")
              
seu <- subset(seu, subset = nFeature_RNA > 250 & nCount_RNA > 500 & percent_mito < 20)  
 geom_vline(xintercept=250, col="black") +
 geom_hline(yintercept=500, col="black")

# Data Normalization
seu <- NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify the most variable genes
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
top10 <- head(VariableFeatures(seu), 10)
write.csv(top10, file = "top10.csv")
plot1 <- VariableFeaturePlot(seu)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 
                 
all.genes <- rownames(seu)

# Scale the data                 
seu <- ScaleData(seu, features = all.genes)

seu <- RunPCA(seu, features = VariableFeatures(object = seu))

# Determine the Dimensionality of Data
                 
ElbowPlot(seu)
VizDimLoadings(seu, dims = 1:2, reduction = "pca")

# Visualize PCA results in a few different ways
DimPlot(seu, reduction = "pca")
DimHeatmap(seu, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(seu, dims = 1:15, cells = 500, balanced = TRUE)
Idents(seu) = "seurat_clusters"

#Determine the dimensionality of the dataset
seu <- JackStraw(seu, num.replicate = 100)
seu <- ScoreJackStraw(seu, dims = 1:20)
JackStrawPlot(seu, dims = 1:20)


seu <- FindNeighbors(seu, dims = 1:20)
seu <- FindClusters(seu, resolution = 1)
ahead(Idents(seu), 5)

# Run UMAP
seu <- RunUMAP(seu, dims = 1:25)     

# Visualize UMAP
DimPlot(seu_diet_merged, reduction = "umap", label = TRUE) + theme(aspect.ratio = 1)

# Visualize UMAPs with different identities
UMAP_sample <- DimPlot(seu_diet2, reduction = "umap", group.by = "orig.ident") + theme(aspect.ratio = 1)
UMAP_sample

UMAP_cohort <- DimPlot(seu_diet2, reduction = "umap", group.by = "cohort") + theme(aspect.ratio = 1)
UMAP_cohort

UMAP_status <- DimPlot(seu, reduction = "umap", group.by = "status") + theme(aspect.ratio = 1)
UMAP_status

# Split the UMAPs
UMAP2 <- DimPlot(seu, reduction = "umap", group.by = "orig.ident", split.by = "cohort") + theme(aspect.ratio = 1)
UMAP2
UMAP3 <- DimPlot(seu, reduction = "umap", group.by = "orig.ident", split.by = "status") + theme(aspect.ratio = 1)
UMAP3


# Add variables to metadata

# Add library type either as enriched for T cells or not
seu$library_type <- case_when(grepl("CD3", seu$orig.ident) ~ "enriched_CD3",
                              grepl("MNC", seu$orig.ident) ~ "MNC")                 

# Add cohort identities
seu$cohorts <- case_when(grepl("9596|2737|2379|2434|2518|4618|6174|9931|1953|25809", seu$orig.ident) ~ "2-Relapsed",
                         grepl("2446|25802|2645|1972|2220|2621|9185|2599", seu$orig.ident) ~ "1-Non-relapsed",
                         grepl("1677|5641|1732|1811|1195|1347|1285|6244|9355|1013", seu$orig.ident) ~ "3-Early relapsed")


# Add timepoints of samples to metadata 
seu $status <- case_when(grepl("9596|2379|4618|9355|2446|1972|1677|1195|5641", seu$orig.ident) ~ "pre_transplant",
                         grepl("25809|2434|2518|6174|9931|1013|25802|2645|2220|2621|9185|2599|1732|1285|6244", seu$orig.ident) ~ "remission",
                         grepl("2737|1953|1811|1347", seu$orig.ident) ~ "relapse")

# Add patient number
seu$patient_identity <- case_when(grepl("2446|25802|2645", seu$orig.ident) ~ "pt01",
                                       grepl("1972|2220|2621", seu$orig.ident) ~ "pt02",
                                       grepl("9185", seu$orig.ident) ~ "pt03",
                                       grepl("2599", seu$orig.ident) ~ "pt04",
                                       grepl("9596|25809|2737", seu$orig.ident) ~ "pt05",
                                       grepl("2379|2434", seu$orig.ident) ~ "pt06",
                                       grepl("2518", seu$orig.ident) ~ "pt07",
                                       grepl("4618|6174|9931|1953", seu$orig.ident) ~ "pt08",
                                       grepl("1677|1732|1811", seu$orig.ident) ~ "pt09",
                                       grepl("1195|1285|1347", seu$orig.ident) ~ "pt10",
                                       grepl("5641|6244", seu$orig.ident) ~ "pt11",
                                       grepl("9355|1013", seu$orig.ident) ~ "pt12")
# Add individual sample ID                

seu$id <- case_when( grepl("2446", seu$orig.ident) ~ "P01.0pre",
                                        grepl("25802_MNC", seu$orig.ident) ~ "P01.1Rem",
                                        grepl("25802_CD3", seu$orig.ident) ~ "P01.1RemT",
                                        grepl("2645", seu$orig.ident) ~ "P01.2Rem",
                                        grepl("1972_MNC", seu$orig.ident) ~ "P02.0pre",
                                        grepl("1972_CD3", seu$orig.ident) ~ "P02.0preT",
                                        grepl("2220", seu$orig.ident) ~ "P02.1Rem",
                                        grepl("2621_MNC", seu$orig.ident) ~ "P02.2Rem",
                                        grepl("2621_CD3", seu$orig.ident) ~ "P02.2RemT",
                                        grepl("9185_MNC", seu$orig.ident) ~ "P03.1Rem",
                                        grepl("9185_CD3", seu$orig.ident) ~ "P03.1RemT",
                                        grepl("2599_MNC", seu$orig.ident) ~ "P04.1Rem",
                                        grepl("2599_CD3", seu$orig.ident) ~ "P04.1RemT",
                                        grepl("9596_MNC", seu$orig.ident) ~ "P05.0pre",
                                        grepl("9596_CD3", seu$orig.ident) ~ "P05.0preT",
                                        grepl("25809", seu$orig.ident) ~ "P05.1Rem",
                                        grepl("2737_MNC", seu$orig.ident) ~ "P05.Rel",
                                        grepl("2737_CD3", seu$orig.ident) ~ "P05.RelT",
                                        grepl("2379_MNC", seu$orig.ident) ~ "P06.0pre",
                                        grepl("2379_CD3", seu$orig.ident) ~ "P06.0preT",
                                        grepl("2434", seu$orig.ident) ~ "P06.1Rem",
                                        grepl("2518_MNC", seu$orig.ident) ~ "P07.1Rem",
                                        grepl("2518_CD3", seu$orig.ident) ~ "P07.1RemT",
                                        
                                        grepl("4618", seu$orig.ident) ~ "P08.0pre",
                                        grepl("6174_MNC", seu$orig.ident) ~ "P08.1Rem",
                                        grepl("6174_CD3", seu$orig.ident) ~ "P08.1RemT",
                                        grepl("9931_MNC", seu$orig.ident) ~ "P08.2Rem",
                                        grepl("9931_CD3", seu$orig.ident) ~ "P08.2RemT",
                                        grepl("1953_MNC", seu$orig.ident) ~ "P08.0Rel",
                                        grepl("1953_CD3", seu$orig.ident) ~ "P08.0RelT",
                                        
                                        grepl("1677_MNC", seu$orig.ident) ~ "P09.0pre",
                                        grepl("1677_CD3", seu$orig.ident) ~ "P09.0preT",
                                        grepl("1732", seu$orig.ident) ~ "P09.1Rem",
                                        grepl("1811_MNC", seu$orig.ident) ~ "P09.Rel",
                                        grepl("1811_CD3", seu$orig.ident) ~ "P09.RelT",
                                        
                                        grepl("1195", seu$orig.ident) ~ "P10.0pre",
                                        grepl("1285", seu$orig.ident) ~ "P10.1Rem",
                                        grepl("1347_MNC", seu$orig.ident) ~ "P10.Rel",
                                        grepl("1347_CD3", seu$orig.ident) ~ "P10.RelT",
                                        
                                        grepl("5641", seu$orig.ident) ~ "P11.0pre",
                                        grepl("6244_MNC", seu$orig.ident) ~ "P11.1Rem",
                                        grepl("6244_CD3", seu$orig.ident) ~ "P11.1RemT",
                                        
                                        grepl("9355", seu$orig.ident) ~ "P12.0pre",
                                        grepl("1013", seu$orig.ident) ~ "P12.1Rem")


seu_diet$timepoint <- case_when(        grepl("2446", seu_diet_merged$orig.ident) ~ "pre-transplant",
                                               grepl("25802", seu_diet_merged$orig.ident) ~ "3mo",
                                               grepl("2645", seu_diet_merged$orig.ident) ~ "6mo",
                                               grepl("1972", seu_diet_merged$orig.ident) ~ "pre-transplant",
                                               grepl("2220", seu_diet_merged$orig.ident) ~ "5mo",
                                               grepl("2621", seu_diet_merged$orig.ident) ~ "2years",
                                               grepl("9185", seu_diet_merged$orig.ident) ~ "1year",
                                               grepl("2599", seu_diet_merged$orig.ident) ~ "3mo",
                                               
                                               
                                               
                                               grepl("9596", seu_diet_merged$orig.ident) ~ "pre-transplant",
                                               grepl("25809", seu_diet_merged$orig.ident) ~ "3mo",
                                               grepl("2737", seu_diet_merged$orig.ident) ~ "relapse-9mo",
                                               grepl("2379", seu_diet_merged$orig.ident) ~ "pre-transplant",
                                               grepl("2434", seu_diet_merged$orig.ident) ~ "3mo",
                                               grepl("2518", seu_diet_merged$orig.ident) ~ "3mo",
                                               grepl("4618", seu_diet_merged$orig.ident) ~ "pre-transplant",
                                               grepl("6174", seu_diet_merged$orig.ident) ~ "3mo",
                                               grepl("9931", seu_diet_merged$orig.ident) ~ "1year", 
                                               grepl("1953", seu_diet_merged$orig.ident) ~ "relapse_2years",
                                               
                                               
                                               grepl("1677", seu_diet_merged$orig.ident) ~ "pre-transplant",
                                               grepl("1732", seu_diet_merged$orig.ident) ~ "1mo",
                                               grepl("1811", seu_diet_merged$orig.ident) ~ "relapse_3mo",
                                               grepl("1195", seu_diet_merged$orig.ident) ~ "pre-transplant",
                                               grepl("1285", seu_diet_merged$orig.ident) ~ "1.5mo",
                                               grepl("1347", seu_diet_merged$orig.ident) ~ "relapse_3mo", 
                                               grepl("5641", seu_diet_merged$orig.ident) ~ "pre-transplant",
                                               grepl("6244", seu_diet_merged$orig.ident) ~ "1.5mo",
                                               grepl("9355", seu_diet_merged$orig.ident) ~ "pre-transplant",
                                               grepl("1013", seu_diet_merged$orig.ident) ~ "1.5mo")


Tcells$timepoint <- case_when(grepl("2446", Tcells$orig.ident) ~ "pre-transplant",
                                         grepl("25802", Tcells$orig.ident) ~ "postTx_3-6m",
                                         grepl("2645", Tcells$orig.ident) ~ "postTx_3-6m",
                                         grepl("1972", Tcells$orig.ident) ~ "pre-transplant",
                                         grepl("2220", Tcells$orig.ident) ~ "postTx_3-6m",
                                         grepl("2621", Tcells$orig.ident) ~ "rem>6m",
                                         grepl("9185", Tcells$orig.ident) ~ "rem>6m",
                                         grepl("2599", Tcells$orig.ident) ~ "postTx_3-6m",
                                         
                                         
                                         
                                         grepl("9596", Tcells$orig.ident) ~ "pre-transplant",
                                         grepl("25809", Tcells$orig.ident) ~ "postTx_3-6m",
                                         grepl("2737", Tcells$orig.ident) ~ "relapse",
                                         grepl("2379", Tcells$orig.ident) ~ "pre-transplant",
                                         grepl("2434", Tcells$orig.ident) ~ "postTx_3-6m",
                                         grepl("2518", Tcells$orig.ident) ~ "postTx_3-6m",
                                         grepl("4618", Tcells$orig.ident) ~ "pre-transplant",
                                         grepl("6174", Tcells$orig.ident) ~ "postTx_3-6m",
                                         grepl("9931", Tcells$orig.ident) ~ "rem>6m", 
                                         grepl("1953", Tcells$orig.ident) ~ "relapse",
                                         
                                         
                                         grepl("1677", Tcells$orig.ident) ~ "pre-transplant",
                                         grepl("1732", Tcells$orig.ident) ~ "rem<3m",
                                         grepl("1811", Tcells$orig.ident) ~ "relapse",
                                         grepl("1195", Tcells$orig.ident) ~ "pre-transplant",
                                         grepl("1285", Tcells$orig.ident) ~ "rem<3m",
                                         grepl("1347", Tcells$orig.ident) ~ "relapse", 
                                         grepl("5641", Tcells$orig.ident) ~ "pre-transplant",
                                         grepl("6244", Tcells$orig.ident) ~ "rem<3m",
                                         grepl("9355", Tcells$orig.ident) ~ "pre-transplant",
                                         grepl("1013", Tcells$orig.ident) ~ "rem<3m")


Tcells$patient_identity <- case_when(grepl("2446|25802|2645", Tcells$orig.ident) ~ "pt01",
                                     grepl("1972|2220|2621", Tcells$orig.ident) ~ "pt02",
                                     grepl("9185", Tcells$orig.ident) ~ "pt03",
                                     grepl("2599", Tcells$orig.ident) ~ "pt04",
                                     grepl("9596|25809|2737", Tcells$orig.ident) ~ "pt05",
                                     grepl("2379|2434", Tcells$orig.ident) ~ "pt06",
                                     grepl("2518", Tcells$orig.ident) ~ "pt07",
                                     grepl("4618|6174|9931|1953", Tcells$orig.ident) ~ "pt08",
                                     grepl("1677|1732|1811", Tcells$orig.ident) ~ "pt09",
                                     grepl("1195|1285|1347", Tcells$orig.ident) ~ "pt10",
                                     grepl("5641|6244", Tcells$orig.ident) ~ "pt11",
                                     grepl("9355|1013", Tcells$orig.ident) ~ "pt12")
 




                            
# Convert each variable to a factor
seu$timepoint <- as.factor(seu@meta.data$orig.ident)
seu$orig.ident <- as.factor(seu@meta.data$orig.ident)
seu$cohort <- as.factor(seu@meta.data$cohort)
seu$id <- as.factor(seu@meta.data$id)
seu$status <- as.factor(seu@meta.data$status)
seu$library_type <- as.factor(seu@meta.data$library_type)
    

# Save the Seurat object at this point 
saveRDS(seu, file = "~/p53_seu.rds")                 

# You can also save as a diet object using the code below, which allows you to store big datasets in much smaller space
seu <-readRDS("~/p53_seu.rds")
seu_diet <- DietSeurat(seu, dimreducs = names(seu@reductions))
saveRDS(seu_diet, file ="~/seu_diet.rds")


seu$sample <- case_when( grepl("2446", seu$orig.ident) ~ "P01.0pre",
                     grepl("25802_MNC", seu$orig.ident) ~ "P01.1Rem",
                     grepl("25802_CD3", seu$orig.ident) ~ "P01.1Rem",
                     grepl("2645", seu$orig.ident) ~ "P01.2Rem",
                     grepl("1972_MNC", seu$orig.ident) ~ "P02.0pre",
                     grepl("1972_CD3", seu$orig.ident) ~ "P02.0pre",
                     grepl("2220", seu$orig.ident) ~ "P02.1Rem",
                     grepl("2621_MNC", seu$orig.ident) ~ "P02.2Rem",
                     grepl("2621_CD3", seu$orig.ident) ~ "P02.2Rem",
                     grepl("9185_MNC", seu$orig.ident) ~ "P03.1Rem",
                     grepl("9185_CD3", seu$orig.ident) ~ "P03.1Rem",
                     grepl("2599_MNC", seu$orig.ident) ~ "P04.1Rem",
                     grepl("2599_CD3", seu$orig.ident) ~ "P04.1Rem",
                     grepl("9596_MNC", seu$orig.ident) ~ "P05.0pre",
                     grepl("9596_CD3", seu$orig.ident) ~ "P05.0pre",
                     grepl("25809", seu$orig.ident) ~ "P05.1Rem",
                     grepl("2737_MNC", seu$orig.ident) ~ "P05.Rel",
                     grepl("2737_CD3", seu$orig.ident) ~ "P05.Rel",
                     grepl("2379_MNC", seu$orig.ident) ~ "P06.0pre",
                     grepl("2379_CD3", seu$orig.ident) ~ "P06.0pre",
                     grepl("2434", seu$orig.ident) ~ "P06.1Rem",
                     grepl("2518_MNC", seu$orig.ident) ~ "P07.1Rem",
                     grepl("2518_CD3", seu$orig.ident) ~ "P07.1Rem",
                     
                     grepl("4618", seu$orig.ident) ~ "P08.0pre",
                     grepl("6174_MNC", seu$orig.ident) ~ "P08.1Rem",
                     grepl("6174_CD3", seu$orig.ident) ~ "P08.1Rem",
                     grepl("9931_MNC", seu$orig.ident) ~ "P08.2Rem",
                     grepl("9931_CD3", seu$orig.ident) ~ "P08.2Rem",
                     grepl("1953_MNC", seu$orig.ident) ~ "P08.0Rel",
                     grepl("1953_CD3", seu$orig.ident) ~ "P08.0Rel",
                     
                     grepl("1677_MNC", seu$orig.ident) ~ "P09.0pre",
                     grepl("1677_CD3", seu$orig.ident) ~ "P09.0pre",
                     grepl("1732", seu$orig.ident) ~ "P09.1Rem",
                     grepl("1811_MNC", seu$orig.ident) ~ "P09.Rel",
                     grepl("1811_CD3", seu$orig.ident) ~ "P09.Rel",
                     
                     grepl("1195", seu$orig.ident) ~ "P10.0pre",
                     grepl("1285", seu$orig.ident) ~ "P10.1Rem",
                     grepl("1347_MNC", seu$orig.ident) ~ "P10.Rel",
                     grepl("1347_CD3", seu$orig.ident) ~ "P10.Rel",
                     
                     grepl("5641", seu$orig.ident) ~ "P11.0pre",
                     grepl("6244_MNC", seu$orig.ident) ~ "P11.1Rem",
                     grepl("6244_CD3", seu$orig.ident) ~ "P11.1Rem",
                     
                     grepl("9355", seu$orig.ident) ~ "P12.0pre",
                     grepl("1013", seu$orig.ident) ~ "P12.1Rem")
           

                 
