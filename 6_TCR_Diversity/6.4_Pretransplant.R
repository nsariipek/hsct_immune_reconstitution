#comparing the pre transplant samples only

#Nurefsan Sariipek
# date: "Dec 20th, 2023"

library(scRepertoire)
library(Seurat)
library(randomcoloR)
library(RColorBrewer)
library(ggpubr)
library(tidyverse)
library(janitor)
library(ggrepel)
library(cowplot)
library(TSP)
library(ggnewscale)
library(ggrastr)

# Empty environment
rm(list=ls())

# For Nurefsan:
my_wd <- "/Users/dz855/Dropbox (Partners HealthCare)/ImmuneEscapeTP53/"
# For Peter:
my_wd <- "~/DropboxMGB/Projects/ImmuneEscapeTP53/"

P01_0pre <- read.csv(paste0(my_wd, "Single Cell Data/2446_MNC/vdj_t/filtered_contig_annotations.csv"))
P02_0pre <-  read.csv(paste0(my_wd, "Single Cell Data/1972_MNC/vdj_t/filtered_contig_annotations.csv"))
P02_0preT <- read.csv(paste0(my_wd, "Single Cell Data/1972_CD3/vdj_t/filtered_contig_annotations.csv"))
                      
P05_0pre <-  read.csv(paste0(my_wd, "Single Cell Data/9596_MNC/vdj_t/filtered_contig_annotations.csv"))
P05_0preT <-  read.csv(paste0(my_wd, "Single Cell Data/9596_CD3/vdj_t/filtered_contig_annotations.csv"))
P06_0pre <-  read.csv(paste0(my_wd, "Single Cell Data/2379_MNC/vdj_t/filtered_contig_annotations.csv"))
P06_0preT <-read.csv(paste0(my_wd, "Single Cell Data/2379_CD3/vdj_t/filtered_contig_annotations.csv"))
P08_0pre <- read.csv(paste0(my_wd, "Single Cell Data/4618_MNC/vdj_t/filtered_contig_annotations.csv"))

# Make a list 
contig_list <- list(P01_0pre, P02_0pre, P02_0preT, P05_0pre, P05_0preT, P06_0pre, P06_0preT, P08_0pre)
combined <- combineTCR(contig_list,
                       samples = c("P01_0pre", "P02_0pre", "P02_0preT", "P05_0pre", "P05_0preT", "P06_0pre", "P06_0preT", "P08_0pre"))

# Add variables
combined <- addVariable(combined, variable.name = "cohort",
                        variables = c("cohort1","cohort1","cohort1",
                                      "cohort2","cohort2","cohort2","cohort2","cohort2"))


combined <- addVariable(combined, variable.name = "ptnumber",
                        variables = c("P01","P02","P02","P05","P05","P06","P06","P08"))

# Visualize diversity metrics
clonalDiversity(combined,
                cloneCall = "strict",
                group.by = "cohort",
                x.axis = "sample",
                n.boots = 100) +
  theme(aspect.ratio = 1,
        axis.text.x = element_text(angle = 45, vjust= 1, hjust = 1, size = 10, color = "black"),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 15, color = "black"))


# Load the Seurat object subsetted for T cells
Tcells <- readRDS(paste0(my_wd, "AnalysisNurefsan/TCR data/RDS/Tcellsubset.rds"))

tabyl(Tcells$id)
  
# Keep only annotated T cell clusters (remove NK cells)
Tcells <- subset(x = Tcells, subset = seurat_clusters %in% c(0,1,2,3,4,5,6,7,9,10,11,12,14)) 
# select only pre_transplant samples
Tcells <- subset(x = Tcells, subset = status == "pre_transplant")
# Check the metadata 
#View(Tcells@meta.data)
# Check the T cell types
table(Idents(Tcells))

# Wrangle the metadata to make it compatible with the TCR metadata (combined)
Tcells_meta <- Tcells@meta.data
Tcells_meta$renamedcells <- gsub("-\\d+_\\d+","", rownames(Tcells@meta.data))
Tcells_meta <- Tcells_meta %>% mutate(fullbc = paste0(Tcells_meta$id, "_", renamedcells, "-1"))
Tcells@meta.data <- Tcells_meta

# Find duplicated barcodes and drop from Seurat object
dup_cells <- Tcells_meta$fullbc[duplicated(Tcells_meta$fullbc)]
UniqueBCs <- Tcells_meta[! Tcells_meta$fullbc %in% dup_cells,]$fullbc 
Tcells <- subset(x = Tcells, subset = fullbc %in% UniqueBCs)
Tcells <- RenameCells(Tcells, new.names = UniqueBCs)


# Combine TCR and sc-RNAseq data
Tcells_combined <- combineExpression(combined, Tcells, 
                                     cloneCall = "strict",
                                     proportion = TRUE,
                                     filterNA = T,
                                     cloneSize = c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=500))


#calculate the frequency, setting group by to sample which is combined samples(T-cell enriched and MNC)
clonalDiversity(Tcells_combined,
                cloneCall = "strict",
                group.by = "sample",
                metrics = c("inv.simpson","gini.simpson"),
                skip.boots = TRUE,
                exportTable = T)

                                           
                                                               