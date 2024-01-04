#Checking all of the samples
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

#cohort 1
P01_0pre <- read.csv(paste0(my_wd, "Single Cell Data/2446_MNC/vdj_t/filtered_contig_annotations.csv"))

P01_1Rem <- read.csv(paste0(my_wd, "Single Cell Data/25802_MNC/vdj_t/filtered_contig_annotations.csv"))
P01_1RemT <- read.csv(paste0(my_wd, "Single Cell Data/25802_CD3/vdj_t/filtered_contig_annotations.csv"))
P01_2Rem <- read.csv(paste0(my_wd, "Single Cell Data/2645_MNC/vdj_t/filtered_contig_annotations.csv"))


P02_0pre <-  read.csv(paste0(my_wd, "Single Cell Data/1972_MNC/vdj_t/filtered_contig_annotations.csv"))
P02_0preT <- read.csv(paste0(my_wd, "Single Cell Data/1972_CD3/vdj_t/filtered_contig_annotations.csv"))
P02_1Rem <-  read.csv(paste0(my_wd, "Single Cell Data/2220_MNC/vdj_t/filtered_contig_annotations.csv"))
P02_2Rem <-  read.csv(paste0(my_wd, "Single Cell Data/2621_MNC/vdj_t/filtered_contig_annotations.csv"))
P02_2RemT <- read.csv(paste0(my_wd, "Single Cell Data/2621_CD3/vdj_t/filtered_contig_annotations.csv"))


P03_1Rem <-    read.csv(paste0(my_wd, "Single Cell Data/9185_MNC/vdj_t/filtered_contig_annotations.csv"))
P03_1RemT <-  read.csv(paste0(my_wd, "Single Cell Data/9185_CD3/vdj_t/filtered_contig_annotations.csv"))

P04_1Rem <-    read.csv(paste0(my_wd, "Single Cell Data/2599_MNC/vdj_t/filtered_contig_annotations.csv"))
P04_1RemT <-  read.csv(paste0(my_wd, "Single Cell Data/2599_CD3/vdj_t/filtered_contig_annotations.csv"))


#cohort 2

P05_0pre <-  read.csv(paste0(my_wd, "Single Cell Data/9596_MNC/vdj_t/filtered_contig_annotations.csv"))
P05_0preT <-  read.csv(paste0(my_wd, "Single Cell Data/9596_CD3/vdj_t/filtered_contig_annotations.csv"))

P05_1Rem <-  read.csv(paste0(my_wd, "Single Cell Data/25809_MNC/vdj_t/filtered_contig_annotations.csv"))

P05_Rel <-  read.csv(paste0(my_wd, "Single Cell Data/2737_MNC/vdj_t/filtered_contig_annotations.csv"))
P05_RelT <-  read.csv(paste0(my_wd, "Single Cell Data/2737_CD3/vdj_t/filtered_contig_annotations.csv"))

P06_0pre <-   read.csv(paste0(my_wd, "Single Cell Data/2379_MNC/vdj_t/filtered_contig_annotations.csv"))
P06_0preT <- read.csv(paste0(my_wd, "Single Cell Data/2379_CD3/vdj_t/filtered_contig_annotations.csv"))
P06_1Rem <-   read.csv(paste0(my_wd, "Single Cell Data/2434_MNC/vdj_t/filtered_contig_annotations.csv"))


P07_1Rem <-  read.csv(paste0(my_wd, "Single Cell Data/2518_MNC/vdj_t/filtered_contig_annotations.csv"))
P07_1RemT <-  read.csv(paste0(my_wd, "Single Cell Data/2518_CD3/vdj_t/filtered_contig_annotations.csv"))

P08_0pre <- read.csv(paste0(my_wd, "Single Cell Data/4618_MNC/vdj_t/filtered_contig_annotations.csv"))

P08_1Rem <-  read.csv(paste0(my_wd, "Single Cell Data/6174_MNC/vdj_t/filtered_contig_annotations.csv"))
P08_1RemT <- read.csv(paste0(my_wd, "Single Cell Data/6174_CD3/vdj_t/filtered_contig_annotations.csv"))

P08_2Rem <-  read.csv(paste0(my_wd, "Single Cell Data/9931_MNC/vdj_t/filtered_contig_annotations.csv"))
P08_2RemT <-  read.csv(paste0(my_wd, "Single Cell Data/9931_CD3/vdj_t/filtered_contig_annotations.csv"))

P08_0Rel <-  read.csv(paste0(my_wd, "Single Cell Data/1953_MNC/vdj_t/filtered_contig_annotations.csv"))
P08_0RelT <-  read.csv(paste0(my_wd, "Single Cell Data/1953_CD3/vdj_t/filtered_contig_annotations.csv"))

P09_0pre <- read.csv(paste0(my_wd, "Single Cell Data/1677_MNC/vdj_t/filtered_contig_annotations.csv"))
P09_0preT <- read.csv(paste0(my_wd, "Single Cell Data/1677_CD3/vdj_t/filtered_contig_annotations.csv"))

P09_1Rem <-  read.csv(paste0(my_wd, "Single Cell Data/1732_MNC/vdj_t/filtered_contig_annotations.csv"))

P09_Rel <-  read.csv(paste0(my_wd, "Single Cell Data/1811_MNC/vdj_t/filtered_contig_annotations.csv"))
P09_RelT <-  read.csv(paste0(my_wd, "Single Cell Data/1811_CD3/vdj_t/filtered_contig_annotations.csv"))

P10_0pre <- read.csv(paste0(my_wd, "Single Cell Data/1195_MNC/vdj_t/filtered_contig_annotations.csv"))

P10_1Rem <-  read.csv(paste0(my_wd, "Single Cell Data/1285_MNC/vdj_t/filtered_contig_annotations.csv"))
P10_Rel <- read.csv(paste0(my_wd, "Single Cell Data/1347_MNC/vdj_t/filtered_contig_annotations.csv"))
P10_RelT <- read.csv(paste0(my_wd, "Single Cell Data/1347_CD3/vdj_t/filtered_contig_annotations.csv"))

P11_0pre <-  read.csv(paste0(my_wd, "Single Cell Data/5641_MNC/vdj_t/filtered_contig_annotations.csv"))
P11_1Rem <- read.csv(paste0(my_wd, "Single Cell Data/6244_MNC/vdj_t/filtered_contig_annotations.csv"))
P11_1RemT <- read.csv(paste0(my_wd, "Single Cell Data/6244_CD3/vdj_t/filtered_contig_annotations.csv"))


P12_0pre <-  read.csv(paste0(my_wd, "Single Cell Data/9355_MNC/vdj_t/filtered_contig_annotations.csv"))
P12_1Rem <-  read.csv(paste0(my_wd, "Single Cell Data/1013_MNC/vdj_t/filtered_contig_annotations.csv"))

# Make a list 
contig_list <- list("P01_0pre", "P01_1Rem", "P01_1RemT", "P01_2Rem", "P02_0pre", "P02_0preT", "P02_1Rem", "P02_2Rem", "P02_2RemT", "P03_1Rem", "P03_1RemT", "P04_1Rem","P04_1RemT", "P05_0pre", "P05_0preT", "P05_1Rem", "P05_Rel", "P05_RelT", "P06_0pre", "P06_0preT", "P06_1Rem", "P07_1Rem", "P07_1RemT", "P08_0Rel", "P08_0RelT", "P08_0pre", "P08_1Rem", "P08_1RemT", "P08_2Rem", "P08_2RemT", "P09_0pre", "P09_0preT", "P09_1Rem", "P09_Rel", "P09_RelT", "P10_0pre", "P10_1Rem", "P10_Rel", "P10_RelT","P11_1Rem", "P11_1RemT", "P12_0pre", "P12_1Rem")
 
combined <- combineTCR(contig_list, 
                       samples = c("P01_0pre", "P01_1Rem", "P01_1RemT", "P01_2Rem", "P02_0pre", "P02_0preT", "P02_1Rem", "P02_2Rem", "P02_2RemT", "P03_1Rem", "P03_1RemT", "P04_1Rem","P04_1RemT", "P05_0pre", "P05_0preT", "P05_1Rem", "P05_Rel", "P05_RelT", "P06_0pre", "P06_0preT", "P06_1Rem", "P07_1Rem", "P07_1RemT", "P08_0Rel", "P08_0RelT", "P08_0pre", "P08_1Rem", "P08_1RemT", "P08_2Rem", "P08_2RemT", "P09_0pre", "P09_0preT", "P09_1Rem", "P09_Rel", "P09_RelT", "P10_0pre", "P10_1Rem", "P10_Rel", "P10_RelT", "P11_1Rem", "P11_1RemT", "P12_0pre", "P12_1Rem"))                    
                       
 combined <- addVariable(combined, variable.name = "status", 
                         variables = c("pre", "rem", "rem", "rem", "pre", "pre","rem","rem","rem","rem","rem","rem","rem","pre","pre","rem","rel","rel","pre","pre","rem","rem","rem","pre","rem","rem","rem","rem","rel","rel"))                      
                      
 combined <- addVariable(combined, variable.name = "ptnumber",
                         variables = c("P01","P01-1","P01-1","P01-2","P02","P02","P02-1","P02-2","P02-2","P03","P03","P04","P04","P05","P05","P05-1","P05-2","P05-2", "P06","P06", "P06-1","P07","P07","P08","P08-1","P08-1","P08-2","P08-2","P08-3","P08-3"))

 combined <- addVariable(combined, variable.name = "cohort",
                         variables = c("cohort1","cohort1","cohort1","cohort1","cohort1","cohort1","cohort1","cohort1","cohort1","cohort1","cohort1","cohort1","cohort1","cohort2","cohort2","cohort2","cohort2","cohort2","cohort2","cohort2","cohort2","cohort2","cohort2","cohort2","cohort2","cohort2","cohort2","cohort2","cohort2","cohort2")) 
 
 # Visualize diversity metrics
 clonalDiversity(combined,
                 cloneCall = "strict",
                 group.by = "ptnumber",
                 x.axis = "status",
                 n.boots = 100) +
   theme(aspect.ratio = 1,
         axis.text.x = element_text(angle = 45, vjust= 1, hjust = 1, size = 10, color = "black"),
         axis.title.x = element_blank(),
         axis.text.y = element_text(size = 15),
         axis.title.y = element_text(size = 15, color = "black"))
 
 # Load the Seurat object subsetted for T cells
 Tcells <- readRDS(paste0(my_wd, "AnalysisNurefsan/TCR data/RDS/Tcellsubset.rds"))
 
 # Keep only annotated T cell clusters (remove NK cells)
 Tcells <- subset(x = Tcells, subset = seurat_clusters %in% c(0,1,2,3,4,5,6,7,9,10,11,12,14)) 
 # Check the metadata 
 #View(Tcells@meta.data)
 # Check the T cell types
 table(Idents(Tcells))
 
 tabyl(Tcells@meta.data$id)
 
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
                                      proportion = FALSE,
                                      group.by = "sample",
                                      filterNA = T,
                                      cloneSize = c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=500))

 
 subset1 <- subsetClones(combined, 
                         name = "sample", 
                         variables = c("P18L", "P18B")) 
 
 
 
 
   
                    
                       