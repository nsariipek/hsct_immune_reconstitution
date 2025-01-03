# Nurefsan Sariipek, Updated at 241231
# Analysis of Immune Escape project
# Load the needed libraries
library(tidyverse)
library(Seurat)
library(ggplot2)
library(readxl)
library(data.table)
library(janitor)
library(cloudml)
library(future.apply)
library(Seurat.utils)

# Use multiple cores
plan(multisession) 

# Start with a clean slate
rm(list=ls())

# Parameters to interact with Google bucket, this part only needed for Terra 
project <- Sys.getenv('WORKSPACE_NAMESPACE')
workspace <- Sys.getenv('WORKSPACE_NAME')
bucket <- Sys.getenv('WORKSPACE_BUCKET')
data_dir <- gs_data_dir_local(bucket)

# Load matrices of all samples
Samples <- c("P1013_MNC","P1732_MNC","P1953_MNC","P2434_MNC", "P2599_CD3","P2791_MNC","P6174_CD3","P9931_CD3",
             "P1195_MNC","P1745_MNC","P1964_MNC","P2446_MNC", "P2599_MNC","P2820_MIX","P6174_MNC","P9931_MNC",
             "P1285_MNC", "P1762_MIX","P1972_CD3","P2448_MNC", "P2621_CD3","P2961_MNC","P6244_CD3",
             "P1347_CD3","P1764_MIX","P1972_MNC","P2517_MIX", "P2621_MNC","P2977_MIX","P6244_MNC",
             "P1347_MNC","P1804_MNC","P2220_MNC","P2518_CD3", "P2645_MNC","P2986_MNC","P9185_CD3",
             "P1665_MIX","P1811_CD3","P2332_MNC","P2518_MNC", "P2698_MIX","P2988_MNC","P9185_MNC",
             "P1671_MIX","P1811_MNC","P2379_CD3","P25809_MNC","P2737_CD3","P3000_MIX","P9355_MNC",
             "P1677_CD3","P1817_MIX","P2379_MNC","P2580_CD3", "P2737_MNC","P4618_MNC","P9596_CD3",
             "P1677_MNC","P1953_CD3","P2408_MNC","P2580_MNC", "P2745_MNC","P5641_MNC","P9596_MNC")



matrices_ls <- future_lapply(Samples, function(x) {
  print(x)
  result <- Read10X((data.dir = paste0("gs/fc-3783b423-62ac-4c69-8c2f-98cb0ee4503b/", x, "/sample_filtered_feature_bc_matrix/")))
  gc()
  return(result)
})


# Create Seurat objects with orig.ident metadata in parallel
seu_ls <- future_lapply(seq_along(matrices_ls), function(i) {
  print(Samples[i])
  result <- CreateSeuratObject(matrices_ls[[i]])
  result@meta.data$orig.ident <- Samples[i]  # Add orig.ident directly
  gc()
  return(result)
})

#Turn to normal version
plan(sequential)

#Create the seurat object combining each sample, add.cell.ids will prevent any duplication due to the different run(same cell identifiers)
seu <- merge(seu_ls[[1]], seu_ls[2:length(seu_ls)], add.cell.ids = Samples)

# Add QC Metrics
seu <- PercentageFeatureSet(seu, "^MT-", col.name = "percent_mito")
seu <- PercentageFeatureSet(seu, "^RB[SL]", col.name = "percent_ribo")

# View Seurat metadata 
View(seu@meta.data)

# Count the cells before filtering 
table(seu$orig.ident)

# QC Control Thresholds
seu <- subset(seu, subset = nFeature_RNA > 250 & nCount_RNA > 500 & percent_mito < 20)   

# Count the cells after filtering
table(seu$orig.ident)

# Add variables to metadata
# Add library type 
seu$library_type <- case_when(grepl("CD3", seu$orig.ident) ~ "CD3",
                              grepl("MNC", seu$orig.ident) ~ "MNC",
                              grepl("MIX", seu$orig.ident) ~ "MIX")                 

# Add cohort information
seu$cohorts <- case_when(grepl("9596|2737|2379|2434|2518|4618|6174|9931|1953|25809", seu$orig.ident) ~ "1-Relapsed",
                         grepl("2446|25802|2645|1972|2220|2621|9185|2599", seu$orig.ident) ~ "1-Non-relapsed",
                         grepl("1677|5641|1732|1811|1195|1347|1285|6244|9355|1013", seu$orig.ident) ~ "1-Early relapsed",
                         grepl("1764|1804|1964|2332|2448|2745", seu$orig.ident) ~ "2-Relapsed",
                         grepl("1665|1745|1817|2408|2988|1762|2698|2791|2977|2986|1671|2517|2820|2961|3000", seu$orig.ident) ~ "2-Non-relapsed")


# Add sample status of samples  
seu$sample_status <- case_when(grepl("9596|2379|4618|9355|2446|1972|1677|1195|5641", seu$orig.ident) ~ "pre_transplant",
                               grepl("25809|2434|2518|6174|9931|1013|25802|2645|2220|2621|9185|2599|1732|1285|6244|1764|1804|1964|2332|2448|2745|1665|1745|1817|2408|2988|1762|2698|2791|2977|2986|1671|2517|2820|2961|3000", seu$orig.ident) ~ "remission",
                               grepl("2737|1953|1811|1347", seu$orig.ident) ~ "relapse"
)

# Add patient number
seu$patient_id <- case_when(grepl("2446|25802|2645", seu$orig.ident) ~ "pt01",
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
                            grepl("9355|1013", seu$orig.ident) ~ "pt12",
                            grepl("1665",seu$orig.ident) ~ "pt13",
                            grepl("1745",seu$orig.ident) ~ "pt14",
                            grepl("1817",seu$orig.ident) ~ "pt15",
                            grepl("2408",seu$orig.ident) ~ "pt16",
                            grepl("2988",seu$orig.ident) ~ "pt17",
                            grepl("1762",seu$orig.ident) ~ "pt18",
                            grepl("2698",seu$orig.ident) ~ "pt19",
                            grepl("2791",seu$orig.ident) ~ "pt20",
                            grepl("2977",seu$orig.ident) ~ "pt21",
                            grepl("2986",seu$orig.ident) ~ "pt22",
                            grepl("1671",seu$orig.ident) ~ "pt23",
                            grepl("2517",seu$orig.ident) ~ "pt24",
                            grepl("2820",seu$orig.ident) ~ "pt25",
                            grepl("2961",seu$orig.ident) ~ "pt26",
                            grepl("3000",seu$orig.ident) ~ "pt27",
                            grepl("1764",seu$orig.ident) ~ "pt28",
                            grepl("1804",seu$orig.ident) ~ "pt29",
                            grepl("1964",seu$orig.ident) ~ "pt30",
                            grepl("2332",seu$orig.ident) ~ "pt31",
                            grepl("2448",seu$orig.ident) ~ "pt32",
                            grepl("2745",seu$orig.ident) ~ "pt33" )


# Add the timepoint to show the samples time as months after tx
seu$timepoint <- case_when(grepl("2446", seu$orig.ident) ~ "pre-transplant",
                                grepl("25802", seu$orig.ident) ~ "3",
                                grepl("2645", seu$orig.ident) ~ "6",
                                grepl("1972", seu$orig.ident) ~ "0",
                                grepl("2220", seu$orig.ident) ~ "5",
                                grepl("2621", seu$orig.ident) ~ "24",
                                grepl("9185", seu$orig.ident) ~ "12",
                                grepl("2599", seu$orig.ident) ~ "3",
                                
                                
                                
                                grepl("9596", seu$orig.ident) ~ "0",
                                grepl("25809", seu$orig.ident) ~ "3",
                                grepl("2737", seu$orig.ident) ~ "9",
                                grepl("2379", seu$orig.ident) ~ "0",
                                grepl("2434", seu$orig.ident) ~ "3",
                                grepl("2518", seu$orig.ident) ~ "3",
                                grepl("4618", seu$orig.ident) ~ "0",
                                grepl("6174", seu$orig.ident) ~ "3",
                                grepl("9931", seu$orig.ident) ~ "12", 
                                grepl("1953", seu$orig.ident) ~ "24",
                                
                                
                                grepl("1677", seu$orig.ident) ~ "0",
                                grepl("1732", seu$orig.ident) ~ "1",
                                grepl("1811", seu$orig.ident) ~ "3",
                                grepl("1195", seu$orig.ident) ~ "0",
                                grepl("1285", seu$orig.ident) ~ "1.5",
                                grepl("1347", seu$orig.ident) ~ "3", 
                                grepl("5641", seu$orig.ident) ~ "0",
                                grepl("6244", seu$orig.ident) ~ "1.5",
                                grepl("9355", seu$orig.ident) ~ "0",
                                grepl("1013", seu$orig.ident) ~ "1.5",
                                grepl("1665",seu$orig.ident) ~ "3",
                                grepl("1745",seu$orig.ident) ~ "3",
                                grepl("1817",seu$orig.ident) ~ "3",
                                grepl("2408",seu$orig.ident) ~ "3",
                                grepl("2988",seu$orig.ident) ~ "3",
                                grepl("1762",seu$orig.ident) ~ "3",
                                grepl("2698",seu$orig.ident) ~ "3",
                                grepl("2791",seu$orig.ident) ~ "3",
                                grepl("2977",seu$orig.ident) ~ "3",
                                grepl("2986",seu$orig.ident) ~ "3",
                                grepl("1671",seu$orig.ident) ~ "3",
                                grepl("2517",seu$orig.ident) ~ "3",
                                grepl("2820",seu$orig.ident) ~ "3",
                                grepl("2961",seu$orig.ident) ~ "3",
                                grepl("3000",seu$orig.ident) ~ "3",
                                grepl("1764",seu$orig.ident) ~ "3",
                                grepl("1804",seu$orig.ident) ~ "3",
                                grepl("1964",seu$orig.ident) ~ "3",
                                grepl("2332",seu$orig.ident) ~ "3",
                                grepl("2448",seu$orig.ident) ~ "3",
                                grepl("2745",seu$orig.ident) ~ "3" )




# Add an unique sample identifier to combine MNC and CD3 libraries
seu$sample_id <- case_when( grepl("2446", seu$orig.ident) ~ "P01_Pre",
                            grepl("25802_MNC", seu$orig.ident) ~ "P01_Rem",
                            grepl("25802_CD3", seu$orig.ident) ~ "P01_Rem",
                            grepl("2645", seu$orig.ident) ~ "P01_Rem",
                            grepl("1972_MNC", seu$orig.ident) ~ "P02_Pre",
                            grepl("1972_CD3", seu$orig.ident) ~ "P02_Pre",
                            grepl("2220", seu$orig.ident) ~ "P02_Rem",
                            grepl("2621_MNC", seu$orig.ident) ~ "P02_Rem",
                            grepl("2621_CD3", seu$orig.ident) ~ "P02_Rem",
                            grepl("9185_MNC", seu$orig.ident) ~ "P03_Rem",
                            grepl("9185_CD3", seu$orig.ident) ~ "P03_Rem",
                            grepl("2599_MNC", seu$orig.ident) ~ "P04_Rem",
                            grepl("2599_CD3", seu$orig.ident) ~ "P04_Rem",
                            
                            grepl("9596_MNC", seu$orig.ident) ~ "P05_Pre",
                            grepl("9596_CD3", seu$orig.ident) ~ "P05_Pre",
                            grepl("25809", seu$orig.ident) ~ "P05_Rem",
                            grepl("2737_MNC", seu$orig.ident) ~ "P05_Rel",
                            grepl("2737_CD3", seu$orig.ident) ~ "P05_Rel",
                            grepl("2379_MNC", seu$orig.ident) ~ "P06_Pre",
                            grepl("2379_CD3", seu$orig.ident) ~ "P06_Pre",
                            grepl("2434", seu$orig.ident) ~ "P06_Rem",
                            grepl("2518_MNC", seu$orig.ident) ~ "P07_Rem",
                            grepl("2518_CD3", seu$orig.ident) ~ "P07_Rem",
                            grepl("4618", seu$orig.ident) ~ "P08_Pre",
                            grepl("6174_MNC", seu$orig.ident) ~ "P08_Rem",
                            grepl("6174_CD3", seu$orig.ident) ~ "P08_Rem",
                            grepl("9931_MNC", seu$orig.ident) ~ "P08_Rem",
                            grepl("9931_CD3", seu$orig.ident) ~ "P08_Rem",
                            grepl("1953_MNC", seu$orig.ident) ~ "P08_Rel",
                            grepl("1953_CD3", seu$orig.ident) ~ "P08_Rel",
                            
                            grepl("1677_MNC", seu$orig.ident) ~ "P09_Pre",
                            grepl("1677_CD3", seu$orig.ident) ~ "P09_Pre",
                            grepl("1732", seu$orig.ident) ~ "P09_Rem",
                            grepl("1811_MNC", seu$orig.ident) ~ "P09_Rel",
                            grepl("1811_CD3", seu$orig.ident) ~ "P09_Rel",
                            
                            grepl("1195", seu$orig.ident) ~ "P10_Pre",
                            grepl("1285", seu$orig.ident) ~ "P10_Rem",
                            grepl("1347_MNC", seu$orig.ident) ~ "P10_Rel",
                            grepl("1347_CD3", seu$orig.ident) ~ "P10_Rel",
                            
                            grepl("5641", seu$orig.ident) ~ "P11_Pre",
                            grepl("6244_MNC", seu$orig.ident) ~ "P11_Rem",
                            grepl("6244_CD3", seu$orig.ident) ~ "P11_Rem",
                            
                            grepl("9355", seu$orig.ident) ~ "P12_Pre",
                            grepl("1013", seu$orig.ident) ~ "P12_Rem",
                            grepl("1665",seu$orig.ident) ~ "P13_Rem",
                            grepl("1745",seu$orig.ident) ~ "P14_Rem",
                            grepl("1817",seu$orig.ident) ~ "P15_Rem",
                            grepl("2408",seu$orig.ident) ~ "P16_Rem",
                            grepl("2988",seu$orig.ident) ~ "P17_Rem",
                            grepl("1762",seu$orig.ident) ~ "P18_Rem",
                            grepl("2698",seu$orig.ident) ~ "P19_Rem",
                            grepl("2791",seu$orig.ident) ~ "P20_Rem",
                            grepl("2977",seu$orig.ident) ~ "P21_Rem",
                            grepl("2986",seu$orig.ident) ~ "P22_Rem",
                            grepl("1671",seu$orig.ident) ~ "P23_Rem",
                            grepl("2517",seu$orig.ident) ~ "P24_Rem",
                            grepl("2820",seu$orig.ident) ~ "P25_Rem",
                            grepl("2961",seu$orig.ident) ~ "P26_Rem",
                            grepl("3000",seu$orig.ident) ~ "P27_Rem",
                            grepl("1764",seu$orig.ident) ~ "P28_Rem",
                            grepl("1804",seu$orig.ident) ~ "P29_Rem",
                            grepl("1964",seu$orig.ident) ~ "P30_Rem",
                            grepl("2332",seu$orig.ident) ~ "P31_Rem",
                            grepl("2448",seu$orig.ident) ~ "P32_Rem",
                            grepl("2745",seu$orig.ident) ~ "P33_Rem")



# Convert each variable to a factor
seu$orig.ident <- as.factor(seu@meta.data$orig.ident)
seu$cohort <- as.factor(seu@meta.data$cohort)
seu$timepoint <- as.factor(seu@meta.data$timepoint)
seu$patient_id <- as.factor(seu@meta.data$patient_id)
seu$sample_status <- as.factor(seu@meta.data$sample_status)
seu$library_type <- as.factor(seu@meta.data$library_type)
seu$sample_id <- as.factor(seu@meta.data$sample_id)


# Visualize QC metrics as a violin plot
feats <- c("nFeature_RNA", "nCount_RNA", "percent_mito", "percent_ribo")

VlnPlot(seu, group.by = "patient_id", features = feats, pt.size = 0.1, ncol = 3) + NoLegend()

FeatureScatter(seu, "nCount_RNA", "nFeature_RNA", group.by = "sample_id", pt.size = 0.5) + theme(aspect.ratio = 1) +
  geom_vline(xintercept=250, col="black") +
  geom_hline(yintercept=500, col="black")

gc()

# Data Normalization
seu <- NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify the st variable genes
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
top10 <- head(VariableFeatures(seu), 10)
write.csv(top10, file = "top10.csv")
plot1 <- VariableFeaturePlot(seu)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1
plot2

gc()

# Scale the data  
all.genes <- rownames(seu)
seu <- ScaleData(seu, features = all.genes)

# Perform linear dimensional reduction
seu <- RunPCA(seu, features = VariableFeatures(object = seu))

# Visualize PCA results in a few different ways
DimPlot(seu, reduction = "pca") #this turned empty?
DimHeatmap(seu, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(seu, dims = 1:15, cells = 500, balanced = TRUE)

# Determine the ‘dimensionality’ of the dataset
ElbowPlot(seu)

#I don't know why this is here
Idents(seu) = "seurat_clusters"

# FInd Neighbors and Cluster cells
seu <- FindNeighbors(seu, dims = 1:15)
seu <- FindClusters(seu, resolution = 1)
head(Idents(seu), 5)

# Run UMAP
seu <- RunUMAP(seu, dims = 1:30)     

# Visualize UMAP
DimPlot(seu, reduction = "umap", label = TRUE) + theme(aspect.ratio = 1)

# Visualize UMAPs with different identities
UMAP_sample <- DimPlot(seu, reduction = "umap", group.by = "orig.ident") + theme(aspect.ratio = 1)
UMAP_sample

UMAP_cohort <- DimPlot(seu, reduction = "umap", group.by = "cohort") + theme(aspect.ratio = 1)
UMAP_cohort

UMAP_status <- DimPlot(seu, reduction = "umap", group.by = "sample_status") + theme(aspect.ratio = 1)
UMAP_status

# Split the UMAPs
UMAP2 <- DimPlot(seu, reduction = "umap", group.by = "orig.ident", split.by = "cohort") + theme(aspect.ratio = 1)
UMAP2
UMAP3 <- DimPlot(seu, reduction = "umap", group.by = "orig.ident", split.by = "status") + theme(aspect.ratio = 1)
UMAP3
#To clear the memory
gc()
# Save the Seurat object at this point 
saveRDS(seu, file = "~/seu.rds")  
#next time you can just load your RDS assigning to a object 
seu <-readRDS("~/seu.rds")

# You can also save as a diet object using the code below, which allows you to store big datasets in much smaller space
seu <-readRDS("~/seu.rds")
seu_diet <- DietSeurat(seu, dimreducs = names(seu@reductions))
saveRDS(seu_diet, file ="~/seu_diet.rds")
# Another method to reduce the size of the obj
seu_filtered <- removeLayersByPattern(seu, pattern = "scale.data", perl = TRUE)

#Run Clusters to do cell annotation(this is the most important step and it takes some time+memory)
seu_markers <- FindAllMarkers(seu, min.pct = .3, logfc.threshold = .3)

#Save as tibble for next time
seu_markers_tib <- as_tibble(seu_markers)
write.csv(seu_markers_tib, file = "~/seu_markers_tib.csv")








