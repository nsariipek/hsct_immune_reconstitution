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

# Cohort 1
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


# Cohort 2

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

P10_Rel <- read.csv(paste0(my_wd, "Single Cell Data/1347_MNC/vdj_t/filtered_contig_annotations.csv"))
P10_RelT <- read.csv(paste0(my_wd, "Single Cell Data/1347_CD3/vdj_t/filtered_contig_annotations.csv"))

P11_0pre <-  read.csv(paste0(my_wd, "Single Cell Data/5641_MNC/vdj_t/filtered_contig_annotations.csv"))
P11_1Rem <- read.csv(paste0(my_wd, "Single Cell Data/6244_MNC/vdj_t/filtered_contig_annotations.csv"))
P11_1RemT <- read.csv(paste0(my_wd, "Single Cell Data/6244_CD3/vdj_t/filtered_contig_annotations.csv"))


P12_0pre <-  read.csv(paste0(my_wd, "Single Cell Data/9355_MNC/vdj_t/filtered_contig_annotations.csv"))


# Make a list 
contig_list <- list(P01_0pre, P01_1Rem, P01_1RemT, P01_2Rem, P02_0pre, P02_0preT, P02_1Rem, P02_2Rem, P02_2RemT, P03_1Rem, P03_1RemT, P04_1Rem, P04_1RemT, P05_0pre, P05_0preT, P05_1Rem, P05_Rel, P05_RelT, P06_0pre, P06_0preT, P06_1Rem, P07_1Rem, P07_1RemT, P08_0Rel, P08_0RelT, P08_1Rem, P08_1RemT, P08_2Rem, P08_2RemT, P09_0pre, P09_0preT, P09_1Rem, P09_Rel, P09_RelT, P10_Rel, P10_RelT, P11_0pre, P11_1Rem, P11_1RemT, P12_0pre)
 
combined <- combineTCR(contig_list, 
                       samples = c("P01_0pre", "P01_1Rem", "P01_1RemT", "P01_2Rem", "P02_0pre", "P02_0preT", "P02_1Rem", "P02_2Rem", "P02_2RemT", "P03_1Rem", "P03_1RemT", "P04_1Rem","P04_1RemT", "P05_0pre", "P05_0preT", "P05_1Rem", "P05_Rel", "P05_RelT", "P06_0pre", "P06_0preT", "P06_1Rem", "P07_1Rem", "P07_1RemT", "P08_0Rel", "P08_0RelT", "P08_1Rem", "P08_1RemT", "P08_2Rem", "P08_2RemT", "P09_0pre", "P09_0preT", "P09_1Rem", "P09_Rel", "P09_RelT", "P10_Rel", "P10_RelT"," P11_0pre", "P11_1Rem", "P11_1RemT", "P12_0pre"))                    
                       
 combined <- addVariable(combined, variable.name = "status", 
                         variables = c("pre", "rem", "rem", "rem", "pre", "pre","rem","rem","rem","rem","rem","rem","rem","pre","pre","rem","rel","rel","pre","pre","rem","rem","rem","rel", "rel", "rem","rem","rem","rem","pre","pre","rem","rel","rel","rel","rel","pre","rem","rem","pre"))                     
                      
 combined <- addVariable(combined, variable.name = "ptnumber",
                         variables = c("P01-0","P01-1","P01-1","P01-2","P02-0","P02-0","P02-1","P02-2","P02-2","P03-1","P03-1","P04-1","P04-1","P05-0","P05-0","P05-1","P05-2","P05-2", "P06-0","P06-0", "P06-1","P07-1","P07-1","P08-3","P08-3","P08-1","P08-1","P08-2","P08-2", "P09-0","P09-0","P09-1","P09-2","P09-2","P10-2","P10-2","P11-0","P11-1","P11-1","P12-0"))

 combined <- addVariable(combined, variable.name = "cohort",
                         variables = c("cohort1","cohort1","cohort1","cohort1","cohort1","cohort1","cohort1","cohort1","cohort1","cohort1","cohort1","cohort1","cohort1","cohort2","cohort2","cohort2","cohort2","cohort2","cohort2","cohort2","cohort2","cohort2","cohort2","cohort2","cohort2","cohort2","cohort2","cohort2","cohort2","cohort3","cohort3","cohort3","cohort3","cohort3","cohort3","cohort3","cohort3","cohort3","cohort3","cohort3")) 
 
 # Visualize diversity metrics
 clonalDiversity(combined,
                 cloneCall = "strict",
                 group.by = "cohort",
                 x.axis = "status",
                 n.boots = 100) +
   theme(aspect.ratio = 1,
         axis.text.x = element_text(angle = 45, vjust= 1, hjust = 1, size = 10, color = "black"),
         axis.title.x = element_blank(),
         axis.text.y = element_text(size = 15),
         axis.title.y = element_text(size = 15, color = "black"))
 
 # Load the Seurat object subsetted for T cells
Tcells <- readRDS(paste0(my_wd, "AnalysisNurefsan/TCR data/RDS/Tcellsfinal.rds"))

 
 # Keep only annotated T cell clusters (remove NK cells)
 Tcells <- subset(x = Tcells, subset = seurat_clusters %in% c(0,1,2,3,4,5,6,7,9,10,11,12,14)) 
 # Check the metadata 
 # View(Tcells@meta.data)
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
                                      proportion = FALSE,
                                      group.by = "ptnumber",
                                      filterNA = T,
                                      cloneSize = c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=500))

 # Adding Souporcell information
 # Load the metadata that contains souporcell information
df1 <- read_csv(paste0(my_wd, "/AnalysisNurefsan/Souporcell/output/cohort1-2_souporcell.csv"))
df3 <- read_csv(paste0(my_wd, "/AnalysisNurefsan/Souporcell/output/cohort3_souporcell.csv"))
combined_df <- bind_rows(df1,df3)
 
 # Wrangle the metadata to 
 combined_df$cell = gsub("_.*","", combined_df$cell)
 combined_df$id = gsub("\\.","_",combined_df$id )
 combined_df$cell= paste0(combined_df$id, "_",combined_df$cell)
 
 
 # Left join the 2 metadata
 Tcells_combined_tib <- as_tibble(Tcells_combined@meta.data, rownames = "cell")
 newdf <- Tcells_combined_tib %>% 
   left_join(combined_df, by ="cell") %>% 
   drop_na()
 
 # Add assignment calls to the Seurat metadata
 Tcells_combined <- AddMetaData(Tcells_combined, data.frame(select(newdf, cell, assignment), row.names = "cell"))
 # Check the numbers
 Tcells_combined$assignment %>% tabyl 
 Tcells_combined_tib <- as_tibble(Tcells_combined@meta.data, rownames = "cell")
 
 tabyl(Tcells_combined_tib, patient_identity, assignment)
 
 
 # Save this RDS file which contains all the TCR data+ souporcell information on the cells 
  saveRDS(Tcells_combined, paste0(my_wd,"AnalysisNurefsan/RDS files/Tcells_combined_all.RDS"))
 

 # Subset the dataset
 pt1 <- subset(x= Tcells_combined, subset = sample %in% c("P01.0pre", "P01.1Rem", "P01.2Rem"))
 pt2 <- subset(x= Tcells_combined, subset = sample %in% c("P02.0pre", "P02.1Rem", "P02.2Rem"))
    
 pt8 <- subset(x= Tcells_combined, subset = sample %in% c("P08.0Rel", "P08.0pre", "P08.1Rem", "P08.2Rem")) 
pt10 <- subset(x= Tcells_combined, subset = sample %in% c("P10.Rel", "P10.0pre", "P10.1Rem"))

Cohort1 <- subset(x= Tcells_combined, subset = cohort == "cohort1") 
 
Cohort2 <- subset(x= Tcells_combined, subset = cohort == "cohort2") 

Cohort3 <- subset(x= Tcells_combined, subset = cohort == "cohort3") 

Pre_tx <- subset(x= Tcells_combined, subset = status == "pre_transplant") 

Rem <- subset(x= Tcells_combined, subset = status == "remission") 
  
Rel <-   subset(x= Tcells_combined, subset = status == "relapse") 

cohort12 <- subset(x= Tcells_combined, subset = cohort %in% c("cohort1", "cohort2"))
  
cohort12pre <- subset(x= cohort12, subset = status == "pre_transplant")
cohort12rem <- subset(x= cohort12, subset = status == "remission")
x <-  subset(x= Tcells_combined, subset =id %in% c("P01_1Rem", "P01_1RemT", "P02_1Rem", "P04_1Rem", "P04_1RemT", "P05_1Rem","P06_1Rem", "P07_1Rem", "P07_1RemT", "P08_1Rem", "P08_1RemT"))

 # Calculate the frequency, setting group by to sample which is combined samples(T-cell enriched and MNC)
 clonalDiversity(x,
                 cloneCall = "strict",
                 group.by = "sample",
                 metrics = c("inv.simpson","gini.simpson"),
                 n.boots = 1000,
                 exportTable = TRUE)
 
 
 

 
 
 
   
                    
                       