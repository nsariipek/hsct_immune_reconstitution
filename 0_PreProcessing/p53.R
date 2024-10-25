#####beginning####
#Nurefsan Sariipek, 230405
#Analysis of data of TP53 project
#Load libraries
library(tidyverse)
library(Seurat)
library(ggplot2)
library(harmony)
library(randomcoloR)
library(readxl)
library(data.table)
library(ggforce)
library(cloudml)
library("RColorBrewer")
library(ggpubr)
library(fgsea)
library(msigdbr)
library(BiocManager)
library("janitor")
library(purrr)
library(fossil)
library(ggrepel)
library(ggalluvial)


# Start with a clean slate
rm(list=ls())

# Parameters to interact with Google bucket
project <- Sys.getenv('WORKSPACE_NAMESPACE')
workspace <- Sys.getenv('WORKSPACE_NAME')
bucket <- Sys.getenv('WORKSPACE_BUCKET')
data_dir <- gs_data_dir_local(bucket)

#Load the saved seurat objects
seu_diet <- readRDS("~/seu_diet.rds")
tc_diet<- readRDS("~/Tcellsubset_diet.rds")
seu_diet_merged <- readRDS("~/seu_diet_merged.rds")


tcremonly <-  subset(x= tc_diet, subset =cohort %in% c("cohort1","cohort2"))
tcremonly <- subset(x= tcremonly, subset =status %in% c("remission"))
tcremonly2_6mo <- subset(x=tcremonly, subset = id %in% c("P01.1Rem", "P01.1RemT", "P01.2Rem", "P02.1Rem", "P04.1Rem", "P04.1RemT", "P05.1Rem", "P06.1Rem", "P07.1Rem", "P07.1RemT", "P08.1Rem", "P08.1RemT"))

#subset only 0-6 months post tx samples from cohorty 1-2
ex <- readRDS("~/cohort1_2_rem0-6.rds")

ex <- subset(x= seu_diet_merged, subset =cohort %in% c("cohort1","cohort2"))
ex <- subset(x= ex, subset =status %in% c("remission"))
ex <- subset(x=ex, subset = id %in% c("P01.1Rem", "P01.1RemT", "P01.2Rem", "P02.1Rem", "P04.1Rem", "P04.1RemT", "P05.1Rem", "P06.1Rem", "P07.1Rem", "P07.1RemT", "P08.1Rem", "P08.1RemT"))

#remove the . for the pseudobulk analysis
ex$id =gsub("\\.", "", ex$id)

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

# Add QC Metrics
seu <- PercentageFeatureSet(seu, "^MT-", col.name = "percent_mito")
seu <- PercentageFeatureSet(seu, "^RB[SL]", col.name = "percent_ribo")
#seu <- PercentageFeatureSet(seu, "^HB[^(P)]", col.name = "percent_hb") 

# Count the cells before filtering 
table(seu_diet_merged$id)

# QC Control Thresholds
seu <- subset(seu, subset = nFeature_RNA > 250 & nCount_RNA > 500 & percent_mito < 20)   

#Count the cells after filtering
table(seu_diet$library_type)

#save as diet object to save place 
seu <-readRDS("~/p53_Nurefsan.rds")
seu_diet <- DietSeurat(seu, dimreducs = names(seu@reductions))
saveRDS(tc_diet, file = "~/Tcellsubset_diet.rds")

seu_diet <- readRDS("~/seu_diet.rds")

saveRDS(seu_diet_merged, file= "~/seu_diet_merged.rds")

#Add variables to metadata

#Add a column with cohort identities
seu_diet $cohorts <- case_when(grepl("9596|2737|2379|2434|2518|4618|6174|9931|1953|25809", seu_diet$orig.ident) ~ "2-Relapse",
                         grepl("2446|25802|2645|1972|2220|2621|9185|2599", seu_diet$orig.ident) ~ "1-Remission",
                         grepl("1677|5641|1732|1811|1195|1347|1285|6244|9355|1013", seu_diet$orig.ident) ~ "3-Early relapse")


#Add timepoints of samples as a column to metadata column 
seu $status <- case_when(grepl("9596|2379|4618|9355|2446|1972|1677|1195|5641", seu$orig.ident) ~ "pre_transplant",
                         grepl("25809|2434|2518|6174|9931|1013|25802|2645|2220|2621|9185|2599|1732|1285|6244", seu$orig.ident) ~ "remission",
                         grepl("2737|1953|1811|1347", seu$orig.ident) ~ "relapse")

#Add patient identities
tc_diet$patient_identity <- case_when(grepl("2446|25802|2645", tc_diet$orig.ident) ~ "pt01",
                                       grepl("1972|2220|2621", tc_diet$orig.ident) ~ "pt02",
                                       grepl("9185", tc_diet$orig.ident) ~ "pt03",
                                       grepl("2599", tc_diet$orig.ident) ~ "pt04",
                                       grepl("9596|25809|2737", tc_diet$orig.ident) ~ "pt05",
                                       grepl("2379|2434", tc_diet$orig.ident) ~ "pt06",
                                       grepl("2518", tc_diet$orig.ident) ~ "pt07",
                                       grepl("4618|6174|9931|1953", tc_diet$orig.ident) ~ "pt08",
                                       grepl("1677|1732|1811", tc_diet$orig.ident) ~ "pt09",
                                       grepl("1195|1285|1347", tc_diet$orig.ident) ~ "pt10",
                                       grepl("5641|6244", tc_diet$orig.ident) ~ "pt11",
                                       grepl("9355|1013", tc_diet$orig.ident) ~ "pt12")

seu_diet$id <- case_when(               grepl("2446", seu_diet$orig.ident) ~ "P01.0pre",
                                        grepl("25802_MNC", seu_diet$orig.ident) ~ "P01.1Rem",
                                        grepl("25802_CD3", seu_diet$orig.ident) ~ "P01.1RemT",
                                        grepl("2645", seu_diet$orig.ident) ~ "P01.2Rem",
                                        grepl("1972_MNC", seu_diet$orig.ident) ~ "P02.0pre",
                                        grepl("1972_CD3", seu_diet$orig.ident) ~ "P02.0preT",
                                        grepl("2220", seu_diet$orig.ident) ~ "P02.1Rem",
                                        grepl("2621_MNC", seu_diet$orig.ident) ~ "P02.2Rem",
                                        grepl("2621_CD3", seu_diet$orig.ident) ~ "P02.2RemT",
                                        grepl("9185_MNC", seu_diet$orig.ident) ~ "P03.1Rem",
                                        grepl("9185_CD3", seu_diet$orig.ident) ~ "P03.1RemT",
                                        grepl("2599_MNC", seu_diet$orig.ident) ~ "P04.1Rem",
                                        grepl("2599_CD3", seu_diet$orig.ident) ~ "P04.1RemT",
                                        grepl("9596_MNC", seu_diet$orig.ident) ~ "P05.0pre",
                                        grepl("9596_CD3", seu_diet$orig.ident) ~ "P05.0preT",
                                        grepl("25809", seu_diet$orig.ident) ~ "P05.1Rem",
                                        grepl("2737_MNC", seu_diet$orig.ident) ~ "P05.Rel",
                                        grepl("2737_CD3", seu_diet$orig.ident) ~ "P05.RelT",
                                        grepl("2379_MNC", seu_diet$orig.ident) ~ "P06.0pre",
                                        grepl("2379_CD3", seu_diet$orig.ident) ~ "P06.0preT",
                                        grepl("2434", seu_diet$orig.ident) ~ "P06.1Rem",
                                        grepl("2518_MNC", seu_diet$orig.ident) ~ "P07.1Rem",
                                        grepl("2518_CD3", seu_diet$orig.ident) ~ "P07.1RemT",
                                        
                                        grepl("4618", seu_diet$orig.ident) ~ "P08.0pre",
                                        grepl("6174_MNC", seu_diet$orig.ident) ~ "P08.1Rem",
                                        grepl("6174_CD3", seu_diet$orig.ident) ~ "P08.1RemT",
                                        grepl("9931_MNC", seu_diet$orig.ident) ~ "P08.2Rem",
                                        grepl("9931_CD3", seu_diet$orig.ident) ~ "P08.2RemT",
                                        grepl("1953_MNC", seu_diet$orig.ident) ~ "P08.0Rel",
                                        grepl("1953_CD3", seu_diet$orig.ident) ~ "P08.0RelT",
                                        
                                        grepl("1677_MNC", seu_diet$orig.ident) ~ "P09.0pre",
                                        grepl("1677_CD3", seu_diet$orig.ident) ~ "P09.0preT",
                                        grepl("1732", seu_diet$orig.ident) ~ "P09.1Rem",
                                        grepl("1811_MNC", seu_diet$orig.ident) ~ "P09.Rel",
                                        grepl("1811_CD3", seu_diet$orig.ident) ~ "P09.RelT",
                                        
                                        grepl("1195", seu_diet$orig.ident) ~ "P10.0pre",
                                        grepl("1285", seu_diet$orig.ident) ~ "P10.1Rem",
                                        grepl("1347_MNC", seu_diet$orig.ident) ~ "P10.Rel",
                                        grepl("1347_CD3", seu_diet$orig.ident) ~ "P10.RelT",
                                        
                                        grepl("5641", seu_diet$orig.ident) ~ "P11.0pre",
                                        grepl("6244_MNC", seu_diet$orig.ident) ~ "P11.1Rem",
                                        grepl("6244_CD3", seu_diet$orig.ident) ~ "P11.1RemT",
                                        
                                        grepl("9355", seu_diet$orig.ident) ~ "P12.0pre",
                                        grepl("1013", seu_diet$orig.ident) ~ "P12.1Rem")
#same thing for whole seu object
seu_diet_merged$id <- case_when(               grepl("2446", seu_diet_merged$orig.ident) ~ "P01.0pre",
                                               grepl("25802_MNC", seu_diet_merged$orig.ident) ~ "P01.1Rem",
                                               grepl("25802_CD3", seu_diet_merged$orig.ident) ~ "P01.1RemT",
                                               grepl("2645", seu_diet_merged$orig.ident) ~ "P01.2Rem",
                                               grepl("1972_MNC", seu_diet_merged$orig.ident) ~ "P02.0pre",
                                               grepl("1972_CD3", seu_diet_merged$orig.ident) ~ "P02.0preT",
                                               grepl("2220", seu_diet_merged$orig.ident) ~ "P02.1Rem",
                                               grepl("2621_MNC", seu_diet_merged$orig.ident) ~ "P02.2Rem",
                                               grepl("2621_CD3", seu_diet_merged$orig.ident) ~ "P02.2RemT",
                                               grepl("9185_MNC", seu_diet_merged$orig.ident) ~ "P03.1Rem",
                                               grepl("9185_CD3", seu_diet_merged$orig.ident) ~ "P03.1RemT",
                                               grepl("2599_MNC", seu_diet_merged$orig.ident) ~ "P04.1Rem",
                                               grepl("2599_CD3", seu_diet_merged$orig.ident) ~ "P04.1RemT",
                                               grepl("9596_MNC", seu_diet_merged$orig.ident) ~ "P05.0pre",
                                               grepl("9596_CD3", seu_diet_merged$orig.ident) ~ "P05.0preT",
                                               grepl("25809", seu_diet_merged$orig.ident) ~ "P05.1Rem",
                                               grepl("2737_MNC", seu_diet_merged$orig.ident) ~ "P05.Rel",
                                               grepl("2737_CD3", seu_diet_merged$orig.ident) ~ "P05.RelT",
                                               grepl("2379_MNC", seu_diet_merged$orig.ident) ~ "P06.0pre",
                                               grepl("2379_CD3", seu_diet_merged$orig.ident) ~ "P06.0preT",
                                               grepl("2434", seu_diet_merged$orig.ident) ~ "P06.1Rem",
                                               grepl("2518_MNC", seu_diet_merged$orig.ident) ~ "P07.1Rem",
                                               grepl("2518_CD3", seu_diet_merged$orig.ident) ~ "P07.1RemT",
                                               
                                               grepl("4618", seu_diet_merged$orig.ident) ~ "P08.0pre",
                                               grepl("6174_MNC", seu_diet_merged$orig.ident) ~ "P08.1Rem",
                                               grepl("6174_CD3", seu_diet_merged$orig.ident) ~ "P08.1RemT",
                                               grepl("9931_MNC", seu_diet_merged$orig.ident) ~ "P08.2Rem",
                                               grepl("9931_CD3", seu_diet_merged$orig.ident) ~ "P08.2RemT",
                                               grepl("1953_MNC", seu_diet_merged$orig.ident) ~ "P08.0Rel",
                                               grepl("1953_CD3", seu_diet_merged$orig.ident) ~ "P08.0RelT",
                                               
                                               grepl("1677_MNC", seu_diet_merged$orig.ident) ~ "P09.0pre",
                                               grepl("1677_CD3", seu_diet_merged$orig.ident) ~ "P09.0preT",
                                               grepl("1732", seu_diet_merged$orig.ident) ~ "P09.1Rem",
                                               grepl("1811_MNC", seu_diet_merged$orig.ident) ~ "P09.Rel",
                                               grepl("1811_CD3", seu_diet_merged$orig.ident) ~ "P09.RelT",
                                               
                                               grepl("1195", seu_diet_merged$orig.ident) ~ "P10.0pre",
                                               grepl("1285", seu_diet_merged$orig.ident) ~ "P10.1Rem",
                                               grepl("1347_MNC", seu_diet_merged$orig.ident) ~ "P10.Rel",
                                               grepl("1347_CD3", seu_diet_merged$orig.ident) ~ "P10.RelT",
                                               
                                               grepl("5641", seu_diet_merged$orig.ident) ~ "P11.0pre",
                                               grepl("6244_MNC", seu_diet_merged$orig.ident) ~ "P11.1Rem",
                                               grepl("6244_CD3", seu_diet_merged$orig.ident) ~ "P11.1RemT",
                                               
                                               grepl("9355", seu_diet_merged$orig.ident) ~ "P12.0pre",
                                               grepl("1013", seu_diet_merged$orig.ident) ~ "P12.1Rem")

tc_diet$id <- case_when(               grepl("2446", tc_diet$orig.ident) ~ "P01.0pre",
                                       grepl("25802_MNC", tc_diet$orig.ident) ~ "P01.1Rem",
                                       grepl("25802_CD3", tc_diet$orig.ident) ~ "P01.1RemT",
                                       grepl("2645", tc_diet$orig.ident) ~ "P01.2Rem",
                                       grepl("1972_MNC", tc_diet$orig.ident) ~ "P02.0pre",
                                       grepl("1972_CD3", tc_diet$orig.ident) ~ "P02.0preT",
                                       grepl("2220", tc_diet$orig.ident) ~ "P02.1Rem",
                                       grepl("2621_MNC", tc_diet$orig.ident) ~ "P02.2Rem",
                                       grepl("2621_CD3", tc_diet$orig.ident) ~ "P02.2RemT",
                                       grepl("9185_MNC", tc_diet$orig.ident) ~ "P03.1Rem",
                                       grepl("9185_CD3", tc_diet$orig.ident) ~ "P03.1RemT",
                                       grepl("2599_MNC", tc_diet$orig.ident) ~ "P04.1Rem",
                                       grepl("2599_CD3", tc_diet$orig.ident) ~ "P04.1RemT",
                                       grepl("9596_MNC", tc_diet$orig.ident) ~ "P05.0pre",
                                       grepl("9596_CD3", tc_diet$orig.ident) ~ "P05.0preT",
                                       grepl("25809", tc_diet$orig.ident) ~ "P05.1Rem",
                                       grepl("2737_MNC", tc_diet$orig.ident) ~ "P05.Rel",
                                       grepl("2737_CD3", tc_diet$orig.ident) ~ "P05.RelT",
                                       grepl("2379_MNC", tc_diet$orig.ident) ~ "P06.0pre",
                                       grepl("2379_CD3", tc_diet$orig.ident) ~ "P06.0preT",
                                       grepl("2434", tc_diet$orig.ident) ~ "P06.1Rem",
                                       grepl("2518_MNC", tc_diet$orig.ident) ~ "P07.1Rem",
                                       grepl("2518_CD3", tc_diet$orig.ident) ~ "P07.1RemT",
                                       
                                       grepl("4618", tc_diet$orig.ident) ~ "P08.0pre",
                                       grepl("6174_MNC", tc_diet$orig.ident) ~ "P08.1Rem",
                                       grepl("6174_CD3", tc_diet$orig.ident) ~ "P08.1RemT",
                                       grepl("9931_MNC", tc_diet$orig.ident) ~ "P08.2Rem",
                                       grepl("9931_CD3", tc_diet$orig.ident) ~ "P08.2RemT",
                                       grepl("1953_MNC", tc_diet$orig.ident) ~ "P08.0Rel",
                                       grepl("1953_CD3", tc_diet$orig.ident) ~ "P08.0RelT",
                                       
                                       grepl("1677_MNC", tc_diet$orig.ident) ~ "P09.0pre",
                                       grepl("1677_CD3", tc_diet$orig.ident) ~ "P09.0preT",
                                       grepl("1732", tc_diet$orig.ident) ~ "P09.1Rem",
                                       grepl("1811_MNC", tc_diet$orig.ident) ~ "P09.Rel",
                                       grepl("1811_CD3", tc_diet$orig.ident) ~ "P09.RelT",
                                       
                                       grepl("1195", tc_diet$orig.ident) ~ "P10.0pre",
                                       grepl("1285", tc_diet$orig.ident) ~ "P10.1Rem",
                                       grepl("1347_MNC", tc_diet$orig.ident) ~ "P10.Rel",
                                       grepl("1347_CD3", tc_diet$orig.ident) ~ "P10.RelT",
                                       
                                       grepl("5641", tc_diet$orig.ident) ~ "P11.0pre",
                                       grepl("6244_MNC", tc_diet$orig.ident) ~ "P11.1Rem",
                                       grepl("6244_CD3", tc_diet$orig.ident) ~ "P11.1RemT",
                                       
                                       grepl("9355", tc_diet$orig.ident) ~ "P12.0pre",
                                       grepl("1013", tc_diet$orig.ident) ~ "P12.1Rem")


seu_diet_merged$timepoint <- case_when(        grepl("2446", seu_diet_merged$orig.ident) ~ "pre-transplant",
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




seu_diet$timepoint <- case_when(        grepl("2446", seu_diet$orig.ident) ~ "pre-transplant",
                                        grepl("25802", seu_diet$orig.ident) ~ "3mo",
                                        grepl("2645", seu_diet$orig.ident) ~ "6mo",
                                        grepl("1972", seu_diet$orig.ident) ~ "pre-transplant",
                                        grepl("2220", seu_diet$orig.ident) ~ "5mo",
                                        grepl("2621", seu_diet$orig.ident) ~ "2years",
                                        grepl("9185", seu_diet$orig.ident) ~ "1year",
                                        grepl("2599", seu_diet$orig.ident) ~ "3mo",
                                        
                                        
                                        
                                        grepl("9596", seu_diet$orig.ident) ~ "pre-transplant",
                                        grepl("25809", seu_diet$orig.ident) ~ "3mo",
                                        grepl("2737", seu_diet$orig.ident) ~ "relapse-9mo",
                                        grepl("2379", seu_diet$orig.ident) ~ "pre-transplant",
                                        grepl("2434", seu_diet$orig.ident) ~ "3mo",
                                        grepl("2518", seu_diet$orig.ident) ~ "3mo",
                                        grepl("4618", seu_diet$orig.ident) ~ "pre-transplant",
                                        grepl("6174", seu_diet$orig.ident) ~ "3mo",
                                        grepl("9931", seu_diet$orig.ident) ~ "1year", 
                                        grepl("1953", seu_diet$orig.ident) ~ "relapse_2years",
                                        
                                        
                                        grepl("1677", seu_diet$orig.ident) ~ "pre-transplant",
                                        grepl("1732", seu_diet$orig.ident) ~ "1mo",
                                        grepl("1811", seu_diet$orig.ident) ~ "relapse_3mo",
                                        grepl("1195", seu_diet$orig.ident) ~ "pre-transplant",
                                        grepl("1285", seu_diet$orig.ident) ~ "1.5mo",
                                        grepl("1347", seu_diet$orig.ident) ~ "relapse_3mo", 
                                        grepl("5641", seu_diet$orig.ident) ~ "pre-transplant",
                                        grepl("6244", seu_diet$orig.ident) ~ "1.5mo",
                                        grepl("9355", seu_diet$orig.ident) ~ "pre-transplant",
                                        grepl("1013", seu_diet$orig.ident) ~ "1.5mo")

seu_diet_merged$groups <- case_when(        grepl("2446", seu_diet_merged$orig.ident) ~ "pre-transplant",
                                            grepl("25802", seu_diet_merged$orig.ident) ~ "remission(3-6mo)",
                                            grepl("2645", seu_diet_merged$orig.ident) ~ "remission(3-6mo)",
                                            grepl("1972", seu_diet_merged$orig.ident) ~ "pre-transplant",
                                            grepl("2220", seu_diet_merged$orig.ident) ~ "remission(3-6mo)",
                                            grepl("2621", seu_diet_merged$orig.ident) ~ "remission(>6mo)",
                                            grepl("9185", seu_diet_merged$orig.ident) ~ "remission(>6mo)",
                                            grepl("2599", seu_diet_merged$orig.ident) ~ "remission(3-6mo)",
                                            
                                            
                                            
                                            grepl("9596", seu_diet_merged$orig.ident) ~ "pre-transplant",
                                            grepl("25809", seu_diet_merged$orig.ident) ~ "remission(3-6mo)",
                                            grepl("2737", seu_diet_merged$orig.ident) ~ "relapse",
                                            grepl("2379", seu_diet_merged$orig.ident) ~ "pre-transplant",
                                            grepl("2434", seu_diet_merged$orig.ident) ~ "remission(3-6mo)",
                                            grepl("2518", seu_diet_merged$orig.ident) ~ "remission(3-6mo)",
                                            grepl("4618", seu_diet_merged$orig.ident) ~ "pre-transplant",
                                            grepl("6174", seu_diet_merged$orig.ident) ~ "remission(3-6mo)",
                                            grepl("9931", seu_diet_merged$orig.ident) ~ "remission(>6mo)", 
                                            grepl("1953", seu_diet_merged$orig.ident) ~ "relapse",
                                            
                                            
                                            grepl("1677", seu_diet_merged$orig.ident) ~ "pre-transplant",
                                            grepl("1732", seu_diet_merged$orig.ident) ~ "remission(0-3mo)",
                                            grepl("1811", seu_diet_merged$orig.ident) ~ "relapse",
                                            grepl("1195", seu_diet_merged$orig.ident) ~ "pre-transplant",
                                            grepl("1285", seu_diet_merged$orig.ident) ~ "remission(0-3mo)",
                                            grepl("1347", seu_diet_merged$orig.ident) ~ "relapse", 
                                            grepl("5641", seu_diet_merged$orig.ident) ~ "pre-transplant",
                                            grepl("6244", seu_diet_merged$orig.ident) ~ "remission(0-3mo)",
                                            grepl("9355", seu_diet_merged$orig.ident) ~ "pre-transplant",
                                            grepl("1013", seu_diet_merged$orig.ident) ~ "remission(0-3mo)")






#Add library type either as enriched for T cells or not
seu$library_type <- case_when(grepl("CD3", seu$orig.ident) ~ "enriched_CD3",
                              grepl("MNC", seu$orig.ident) ~ "MNC")
                            
#use this line to convert to factor
seu_diet_merged$groups <- as.factor(seu_diet_merged@meta.data$groups)
seu$cohort <- as.factor(seu@meta.data$cohort)
tc_diet$id <- as.factor(tc_diet@meta.data$id)
seu$library_type <- as.factor(seu@meta.data$library_type)

#check the metadata
View(seu_diet_merged@meta.data)  
View(tc_diet@meta.data)

#use this to visualize the numbers in different cohorts
as_tibble(seu_diet@meta.data) %>% tabyl(id, celltype)
as_tibble(seu_diet_merged@meta.data) %>% tabyl(timepoint)
as_tibble(tc_diet@meta.data) %>% tabyl(celltype)
table(seu_diet$orig.ident)
table(seu$cohort)
table(seu_diet_merged$groups)

# Visualize QC metrics as a violin plot
feats <- c("nFeature_RNA", "nCount_RNA", "percent_mito", "percent_ribo")
VlnPlot(seu_diet, group.by = "patient_identity", features = feats, pt.size = 0.1, ncol = 3) + NoLegend()
FeatureScatter(seu_diet_merged, "nCount_RNA", "nFeature_RNA", group.by = "id", pt.size = 0.5) + theme(aspect.ratio = 1) +
  geom_vline(xintercept=250, col="black") +
  geom_hline(yintercept=500, col="black")

seu <- subset(seu, subset = nFeature_RNA > 250 & nCount_RNA > 500 & percent_mito < 20)  
geom_vline(xintercept=250, col="black") +
geom_hline(yintercept=500, col="black")

#Data Normalization
seu <- NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 10000)

#Identify the most variable genes
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
top10 <- head(VariableFeatures(seu), 10)
write.csv(top10, file = "top10.csv")
plot1 <- VariableFeaturePlot(seu)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 
all.genes <- rownames(seu)

seu <- ScaleData(seu, features = all.genes)

seu <- RunPCA(seu, features = VariableFeatures(object = seu))

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
seu <- FindClusters(seu, resolution = 1)
ahead(Idents(seu), 5)

#Run UMAP
seu <- RunUMAP(seu, dims = 1:25)

#Visualize UMAP
DimPlot(seu_diet_merged, reduction = "umap", label = TRUE) + theme(aspect.ratio = 1)

#Visualize UMAPs with different identities
UMAP_sample <- DimPlot(seu_diet2, reduction = "umap", group.by = "orig.ident") + theme(aspect.ratio = 1)
UMAP_sample

UMAP_cohort <- DimPlot(seu_diet2, reduction = "umap", group.by = "cohort") + theme(aspect.ratio = 1)
UMAP_cohort

UMAP_status <- DimPlot(seu, reduction = "umap", group.by = "status") + theme(aspect.ratio = 1)
UMAP_status

#split the UMAPs
UMAP2 <- DimPlot(seu, reduction = "umap", group.by = "orig.ident", split.by = "cohort") + theme(aspect.ratio = 1)
UMAP2
UMAP3 <- DimPlot(seu, reduction = "umap", group.by = "orig.ident", split.by = "status") + theme(aspect.ratio = 1)
UMAP3

#saveRDS as an object at this point and do not modify it after this part meaning do not re-save it
saveRDS(seu_diet, file = "~/seu_diet.rds")

#Run FindAllMarkers for to see most expressed genes
seu_markers <- FindAllMarkers(seu, min.pct = .3, logfc.threshold = .3)

#Save as tibble and export as a csv for next time and also to export and work on cell type annotations
seu_markers_tib <- as_tibble(seu_markers)
write.csv(seu_markers_tib, file = "~/seu_markers_tib.csv")

#Load marker genes that you have saved before
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
FeaturePlot(seu_diet, features = c("CD19", "MS4A1", "PDCD1LG2", "NT5E", "FCER2", "SDC1",
                              "PAX5", "TCF3", "CD80", "SPIB", "BCL6")) 


FeaturePlot(seu_diet_merged, features = c("CLEC11A", "PRAME", "AZU1", "NREP", "ARMH1", "C1QBP", "TRH"))

###### Rename clusters #####
View(seu_diet@meta.data)

Idents(seu_diet) = "seurat_clusters"

#Rename Clusters
seu.cluster.ids <- c("T cells","T cells","T cells","Monocytes","Blasts","T cells","T cells","Late Erythroids","Monocytes","Mid Erythroids", "B cells","T cells","Pro Monocytes", "Mid Erythroids","Mid Erythroids","Late Erythroids","Blasts","Early Erythroids","Non Classical Monocytes", "B cells", "Monocytes", "Monocytes", "cDC","Early Erythroids","Pre B cells", "HSPCs", "Unidentified", "Plasma Cells", "Blasts", "Blasts", "Pro B cells", "Late Erythroids", "pDC", "Late Erythroids", "Doublets")

names(seu.cluster.ids) <- levels(seu_diet)
seu_diet <- RenameIdents(seu_diet, seu.cluster.ids)
seu_diet@meta.data$celltype = Idents(seu_diet)

#see the levels
levels(seu_diet$celltype)
levels(seu_diet)
View(seu_diet@meta.data)

#Visualize the annotated clusters
mycolors <- distinctColorPalette(k = 34)
pie(rep(1, 34), col = mycolors) 
DimPlot(seu_diet_merged, reduction = "umap", repel = T, group.by = "celltype", label = T) + theme(aspect.ratio = 1)

####Subset of T cells ####

#Load diet seurat object from files 
tc_diet<- readRDS("~/Tcellsubset_diet.rds")

#Load saved Tcell subset seurat object and and turn it to a diet seurat for next times
tc <- readRDS("~/Tcellsubset.RDS")
tc_diet <- DietSeurat(tc, dimreducs = names(tc@reductions))
saveRDS(tc_diet, file = "~/Tcellsubset_diet.rds")
tc_diet<- readRDS("~/Tcellsubset_diet.rds")

Tcell_subset2 <- subset(x = seu, subset = seurat_clusters %in% c(0,1,2,5,6,11))

#Don't forget to run FindVariables again for more accurate analysis
Tcell_subset2 <- FindVariableFeatures(Tcell_subset2, selection.method = "vst", nfeatures = 2000, verbose = FALSE)

#Decide dimensions according to elbow plot below
ElbowPlot(Tcell_subset2)

Tcell_subset2 <- FindNeighbors(Tcell_subset2, dims = 1:20)
Tcell_subset2 <- FindClusters(Tcell_subset2, resolution = 1.1)

#Run UMAP
Tcell_subset2 <- RunUMAP(Tcell_subset2, dims = 1:20) 
DimPlot(Tcell_subset2, reduction = "umap", label = TRUE) + theme(aspect.ratio = 1)

#Save for the next time
saveRDS(Tcell_subset2, file= "~/Tcellsubset.RDS")


#T cell Features to distinguish populations
FeaturePlot(tc_diet, features = c("CD8A", "CD8B", "CD4", "NCAM1","IL10","TGFB","GATA3","TCF7","SELL","CCR7","SELL","TMIGD2","LEF1","CD28","CD27"))

FeaturePlot(tc_diet, features = c("CD8A", "CD8B", "CD4", "NCAM1","IFNG", "CCL3", "PDCD1", "TIGIT", "LAG3", "HAVCR2", "CTLA4"))

#Run Find Markers
Tcell_subset_markers <- FindAllMarkers(Tcell_subset2, min.pct = .3, logfc.threshold = .3)

#Convert to a tibble and export it as a csv
Tcell_subset_markers_tibble <- as_tibble(Tcell_subset_markers)
write.csv(Tcell_subset_markers_tibble, file = "~/Tcell_subset_markers.csv")

#use this library to visualize the numbers in different cohorts
library("janitor")
as_tibble(tc@meta.data) %>% tabyl(status)
as_tibble(tc@meta.data) %>% tabyl(cohort)
as_tibble(tc@meta.data) %>% tabyl(library_type)
table(tc_diet$status)
table(tc_diet$cohort)
table(tc_diet$library_type)


#Visualize the populations
Tcell_UMAP_sample <- DimPlot(tc_diet2, reduction = "umap", group.by = "orig.ident") + theme(aspect.ratio = 1)
Tcell_UMAP_sample
Tcell_UMAP_status <- DimPlot(tc_diet, reduction = "umap", group.by = "status") + theme(aspect.ratio = 1)
Tcell_UMAP_status
Tcell_UMAP_cohort <- DimPlot(tc_diet2, reduction = "umap", group.by = "cohort") + theme(aspect.ratio = 1) 
Tcell_UMAP_cohort

#For annotation purposes upload Kyle's signatures and add as a module score
Kylecells <- read.csv(file = "~/KR_CellTypeSignatures.csv", header = T)

#Add Module Scores 
for (n in names(Kylecells)) {
  print(n)
  #n <- "HSPC"
seu_diet <- AddModuleScore(object = seu_diet, features = Kylecells[n], name = n)
  colnames(seu_diet@meta.data) <- gsub(str_c(n, "1$"), str_c(n, "_Score"), colnames(seu_diet@meta.data))}

#visualize
FeaturePlot(seu_diet, features = c("NK_Score","Monocyte_Score","B_cells_Score","Early_Erythroid_Score", "Mid_Erythroid_Score", "Late_Erythroid_Score", "cDC_Score", "pDC_Score", "Plasma_Cell_Score", "Pre_B_cell_Score", "Cycling_NP_Score", "NP_Score", "ProB_Score", "GMP_Score", "HSPC_Score", "EP_Score", "IMP_Score", "MkP_Score", "MPP_Score", "EBM_Score", "CD4_Naïve_Score", "CD56_dim_NK_Score", "CD8_Term_Eff_Score", "CD8_GZMK_Exh_Score", "CD8_EM_Score", "CD8_Naïve_Score", "NKT_Score", "CD4_CM_Score", "MAIT_Score", "Tregs_Score", "CD56_Bright_NK_Score"))

#For annotation purposes upload Erica's signatures and add as a module score
markergenes <- read.table(file = "~/markerGenes.txt", header = T)

for (n in names(markergenes)) {
  print(n)
  #n <- "HSPC"
  seu_diet <- AddModuleScore(object = seu_diet, features = markergenes[n], name = n)
  colnames(seu_diet@meta.data) <- gsub(str_c(n, "1$"), str_c(n, "_Score"), colnames(seu_diet@meta.data))
}

#Visualize
FeaturePlot(seu_diet, features = c("HSPC_Score",	"EarlyEry_Score",	"LateEry_Score",	
                              "GMP_Score",	"ProMono_Score",	"Mono_Score",	"ncMono_Score",	
                              "cDC_Score",	"pDC_Score",	"ProB_Score",	"PreB_Score",	
                              "B_Score",	"Plasma_Score",	"CD4Naive_Score",	"CD4Memory_Score",
                              "CD8Naive_Score",	"CD8Memory_Score",	"CD8TermExh_Score",
                              "GammaDeltaLike_Score",	"NKT_Score",	"NK_Score"))

#########T cell Annotations ########

#Heatmap the genes in your interest across the clusters#

#Melanoma paper
gene_list = c("CD3E","CD4","CD8A","SELL","CCR7","IL7R","CD28","FAS","CD27","ITGAE","ITGAL","ITGAM","ITGAX","PDCD1","TIGIT","HAVCR2","LAG3","CTLA4","VTCN1","CD244","KLRG1","TNFRSF14","BTLA","CD160","CD38","ENTPD1","NT5E","CD69", "IL2RA","ICOS","TNFRSF4","TNFRSF9","HLA-DRA","CD40LG","GZMA","GZMB","GZMH","GZMK","GZMM","PRF1","NKG7","GNLY","IFNG","FASLG","TNF","IL17A","IL2","LEF1","TCF7","EOMES","TBX21","PRDM1","TOX","GATA3","ID2","ID3","NR4A1","ZNF683","FOXP3","MKI67","TOP2A","TRGV9","TRDV2","KLRB1","KLRC3")

gene_list = c("CD3E","CD4","CD8A","SELL","CCR7","IL7R","CD28","FAS","CD27","ITGAE","ITGAL","ITGAM","ITGAX","PDCD1","TIGIT","HAVCR2","LAG3","CTLA4","VTCN1","CD244","KLRG1","TNFRSF14","BTLA","CD160","CD38","ENTPD1","NT5E","CD69", "IL2RA","ICOS","TNFRSF4","TNFRSF9","HLA-DRA","CD40LG","GZMA","GZMB","GZMH","GZMK","GZMM","PRF1","NKG7","GNLY")


#CD8
gene_list = c("CD8","CD8A","CD8B","CD4","NKG7","GNLY", "CST7", "PRF1", "GZMK", "GZMH","GZMA", "GZMB", "IFNG", "CCL3", "PDCD1", "TIGIT", "LAG3", "HAVCR2", "CTLA4", "TCF7", "LEF1", "SELL", "CD27","CD28","CD57", "S1PR1", "VIM", "GPR183", "CCR7", "IL7R","CCL3", "CCL3L1", "CCL4L2", "CCL4")

#CD4 
gene_list = c("CD4","NKG7","GNLY", "CCL4", "CST7", "PRF1", "GZMK", "GZMH","GZMA", "GZMB", "IFNG", "CCL3", "PDCD1", "TIGIT", "LAG3", "HAVCR2", "CTLA4","FOXP3", "IKZF2", "IL2RA", "CCR10", "TNFRSF4","TIMP1", "USP10", "CCR7", "TCF7", "LEF1", "SELL","IL7R","ANXA1","CD40LG", "RORA")

#CD4
gene_list = c("CCR7","SELL","TMIGD2","LEF1","ANXA1","LGALS1","TIMP1","S100A11","ANXA2","KLBR1","CCL5","GZMA","GZMK")

#CD56
gene_list = c("NCAM", "CD8","CD8A","CD8B", "GZMK", "XCL1", "IL7R", "TCF7", "GPR183", "GZMB", "PRF1", "CX3CR1", "CD7", "FCER1G", "KLRB1","KLRC2", "CD3E", "PATL2", "ZBTB38")

#th1 and 2 
gene_list = c("KLRD1", "IFNGR1","CXCR3","CXCR6","CCR1","CCR5","STAT1","STAT4","TBX21","TNF","LTA","IFNG","IL2","IL12RB1","IL18R1","TNFSF11","HAVCR","CXCR4","BATF","IL4","CCR4","GATA3","IL5","IL13","CCR8","IRF4","AREG","STAT6","HAVCR1","PTGDR2","IL17RB","IL33","IL1R1", "AHR","CSF2","KLRB1","BATF","IL17A", "CCR4", "MAF", "IL17AF","CCR6","NFKBIZ","IL17F","IL21R","IRF4","IL21","IL22","RORA","RORC","STAT3","TBX21","PRF1","GZMB","GZMA")

gene_list = c("IL10","TGFB")

tc_diet_avg <- AverageExpression(tc_diet, return.seurat = TRUE)

tc_diet_avg_cd8 <- AverageExpression(tc_diet %>% subset(seurat_clusters %in% c(1,6)), return.seurat = TRUE)

tc_diet_avg_nk <- AverageExpression(tc_diet %>% subset(seurat_clusters %in% c(8,11,13)), return.seurat = TRUE)

tc_diet_avg_cd4 <- AverageExpression(tc_diet %>% subset(seurat_clusters %in% c(0,3,4,7,10,14)), return.seurat = TRUE)

#Visualize it either as heatmap or dot plot
DoHeatmap(tc_diet_avg, features = gene_list, draw.lines = FALSE)  +  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 5, name = "RdBu"))) + theme_pubr(base_size = 10) + guides(colour = "none")

DotPlot(tc_diet %>% subset(seurat_clusters %in% c(1,2,5,6,9)), features = gene_list, dot.scale = 10) + theme(axis.text.x = element_text(angle=90,vjust = 0.5), axis.text = element_text(size=11), text = element_text(size=11)) + scale_color_distiller(palette = "RdBu")

#Rename T cell Identities
Tcell.cluster.ids <- c("CD4 Memory","CD8 Effector","CD8 Memory","CD4 Naïve","Treg","CD8 Naïve","CD8 Effector","CD4 Naïve","CD56 Dim NK cells","CD8 Terminally Exhausted","CD4 Memory","NK T cells","γδ T lymphocytes","CD56 Bright NK cells", "CD4 Naïve")
names(Tcell.cluster.ids) <- levels(tc_diet)
tc_diet <- RenameIdents(tc_diet, Tcell.cluster.ids)
tc_diet@meta.data$celltype = Idents(tc_diet)

#see the new annotated UMAP
table(seu_diet_merged)
Idents(tc_diet) = "seurat_clusters"
DimPlot(tc_diet, reduction = "umap") + theme(aspect.ratio = 1)

##### Merge cell type annotation#####
meta = seu_diet@meta.data
tnk_meta = tc_diet@meta.data

meta$celltype = as.character(meta$celltype)
tnk_meta$celltype = as.character(tnk_meta$celltype)
meta$cell = rownames(meta)
tnk_meta$cell = rownames(tnk_meta)
meta$celltype[meta$cell %in% tnk_meta$cell] = tnk_meta$celltype
seu_diet@meta.data = meta

#Save the merged seu as a new object
seu_diet_merged <- DietSeurat(seu_diet, dimreducs = names(seu_diet@reductions))
saveRDS(seu_diet_merged, file = "~/seu_diet_merged.rds")
seu_diet_merged <- readRDS("~/seu_diet_merged.rds")

#Visualize the new UMAP
Idents(seu_diet_merged) = "celltype"

DimPlot(seu_diet_merged, reduction = "umap", label = TRUE) + theme(aspect.ratio = 1)

DimPlot(seu_diet_merged, label = F, repel = T, label.size = 4) + theme(text = element_text(size=18)) & guides(color=guide_legend(ncol =1)) + theme(aspect.ratio = 1) 

DimPlot(seu_diet_merged, label = F, repel = T, label.size = 4) & scale_color_manual(values = paletteer_c("grDevices::Set 2", 30)) & theme(text = element_text(size=10)) & guides(color=guide_legend(ncol =1)) + theme(aspect.ratio = 1)

paletteer_c("grDevices::Set 2", 30)

table(Idents(seu_diet_merged))

palette={ "Unidentified": "gray" }

#Data Normalization
seu_diet <- NormalizeData(seu_diet, normalization.method = "LogNormalize", scale.factor = 10000)
gc()

#Identify the most variable genes
seu_diet <- FindVariableFeatures(seu_diet, selection.method = "vst", nfeatures = 2000, verbose = FALSE) 

#Run FindAllMarkers for to see most expressed genes in fully annotated object
seu_markers <- FindAllMarkers(seu_diet, min.pct = .3, logfc.threshold = .3)

#Save as tibble and export as a csv for next time and also to export and work on cell type annotations
seu_markers_tib <- as_tibble(top50)
write.csv(seu_markers_tib, file = "~/seu_markers_tib.csv")

seu_markers %>%
  group_by(cluster) %>%
  top_n(n = 50, wt = avg_log2FC) -> top50

######Proportion graphs#####

## Plot showing cell proportions for whole data set##
MNC_subset@meta.data %>% 
  group_by(celltype, cohorts) %>% 
  summarize(n=n()) %>% 
  ggplot(aes(x=cohorts, y=n, fill=celltype)) + 
  geom_bar(stat = "identity", position = "fill") + 
  theme_pubr(base_size = 15) +
  scale_fill_manual(values=palette) +
  theme(axis.text.x = element_text(angle=90, vjust = 0.5), legend.position = "right") +
  scale_fill_manual(values = c(brewer.pal(12,"Paired"), brewer.pal(12,"Paired"),brewer.pal(12,"Paired"),brewer.pal(12,"Paired")))

#to add a color palette
palette <- distinctColorPalette(k = 34)
pie(rep(1, 34), col = palette) 

## Plot showing cell proportions for T cells ##
tc_diet@meta.data %>% 
  select(celltype, cohort) %>% 
  group_by(celltype, cohort) %>% 
  summarize(n=n()) %>% 
  ggplot(aes(x=cohort, y=n, fill=celltype)) + 
  geom_bar(stat = "identity", position = "fill") + 
  theme_pubr(base_size = 15) +
  scale_fill_manual(values = c(brewer.pal(12,"Paired"), "black", (brewer.pal(12,"Paired")))) +
  theme(axis.text.x = element_text(angle=90, vjust = 0.5), legend.position = "right")

#Cell proportion bar graph per cohorts side to side
#Summarize for plotting purposes
meta_summary <- 
  tc_diet@meta.data %>% group_by(celltype, cohort) %>%
  summarize(n_cells = n())
#add proportion
meta_summary <- meta_summary %>% group_by(cohort) %>% arrange(cohort, -n_cells) %>%
  mutate(n_norm = n_cells/sum(n_cells))

p1 <- meta_summary %>%
  ggplot(aes(x=cohort, y = n_norm, fill=celltype)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_pubr(base_size = 15) +
  ylab("Percentages") +
  ylim(0.0,0.3)+
  scale_fill_manual(values =c("#00FF40CC","#FFA500", "#FF0800CC")) +
  theme(aspect.ratio = 1, axis.text.x = element_text(angle = 45, vjust= 1, hjust = 1, size = 10),  legend.position = "right")

p2 <- meta_summary %>%
  ggplot(aes(x=cohort, y = n_norm, fill=celltype)) +
  geom_bar(stat = "identity", position = "fill") +
  theme_pubr(base_size = 15) +
  ylab("Percentages") +
  scale_fill_manual(values = c(brewer.pal(12,"Paired"), "black", (brewer.pal(12,"Paired")))) +
  theme(axis.text.x = element_text(angle=90, vjust = 0.5), legend.position = "right")
p2

######Subset the libraries ####
#Whole object-  MNC libraries only
MNC_subset <- subset(x = seu_diet, subset = library_type %in% c("MNC"))

#T cell subset- MNC libraries only 
MNC_Tsubset <- subset(x = tc_diet, subset = library_type %in% c("MNC"))

#Subset T enriched libraries, T cells 
CD3_subset <- subset(x = tc_diet, subset = library_type %in% c("enriched_CD3")) 

#Whole data excluding enriched libraries
Pre_transplant <- subset(x = MNC_subset, subset = status %in% c("pre_transplant"))
Remission <- subset(x = MNC_subset, subset = status %in% c("remission"))
Relapse <- subset(x = MNC_subset, subset = status %in% c("relapse"))

#T cell seurat object
T_Pre_transplant <- subset(x = MNC_Tsubset, subset = status %in% c("pre_transplant"))
T_Relapse <- subset(x = MNC_Tsubset, subset = status %in% c("relapse"))
T_Remission <- subset(x = tc_diet, subset = status %in% c("remission"))

#Subset only 2-6 mo remission samples
T_Remission2_6 <- subset(x = tc_diet, subset = id %in% c("P01.1Rem", "P01.1RemT", "P02.1Rem", "P04.1Rem", "P04.1RemT", "P05.1Rem","P06.1Rem", "P07.1Rem", "P07.1RemT", "P08.1Rem", "P08.1RemT"))


view(T_Remission2_6@meta.data)
Tcells_tib <- as_tibble(T_Remission2_6@meta.data) %>% tabyl(celltype, cohort, patient_identity)
write.csv(Tcells_tib, file = "~/remission_2-6mo_Tcells_tib.csv")


#use this library to visualize the numbers in different cohorts
library("janitor")
as_tibble(CD3_subset@meta.data) %>% tabyl(status)
as_tibble(tc_diet@meta.data) %>% tabyl(library_type)
as_tibble(MNC_subset@meta.data) %>% tabyl(library_type)
table(MNC_subset$status)
table(MNC_Tsubset$cohort)
table(tc_diet$library_type)

#Cell proportions plot 
T_Remission2_6@meta.data %>% 
  group_by(celltype, cohorts) %>% 
  summarize(n=n()) %>% 
  ggplot(aes(x=celltype, y=n, fill=cohorts)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  theme_pubr(base_size = 15) +
  theme(axis.text.x = element_text(angle=90, vjust = 0.5), legend.position = "right") +
  scale_fill_manual(values =c("#00FF40CC", "#FF0800CC"))


##alluvial plot

#This code creates a summary metadata that will allow us to see cell proportions over the time. 
#check the metadata
View(seu_diet_merged@meta.data) 

#select the variables you're interested in

t = seu_diet_merged@meta.data %>% rownames_to_column("barcode")%>% 
  filter(library_type=="MNC") %>% 
  filter(cohort == "cohort1")%>% 
  dplyr::select(groups,id,celltype,patient_identity,cohort) 

#see the new metadata 
head(t)

#only select the t cells
t <- subset(t, subset = celltype %in% c("CD4 Memory","CD8 Effector","CD8 Memory","CD4 Naïve","Treg","CD8 Naïve","CD8 Effector","CD4 Naïve","CD56 Dim NK cells","CD8 Terminally Exhausted","CD4 Memory","NK T cells","γδ T lymphocytes","CD56 Bright NK cells", "CD4 Naïve"))

#summarise and calculate the frequencies 
t_sum = t %>% 
  group_by(groups, cohort, celltype) %>%
  summarize(n = n()) %>%
  group_by(groups, cohort) %>%
  mutate(freq = n/sum(n)) %>%
  dplyr::select(-n)

head(t_sum)

library(ggalluvial)
#select the patient 
t_sum1= t_sum %>%  
  mutate(celltype = factor(celltype, level = celltype)) %>%  
  filter(cohort=="cohort2") 

head(t_sum1)

ordered_celltypes = c("CD4 Naïve","CD4 Memory","Treg","CD8 Naïve","CD8 Memory","CD8 Effector","CD8 Terminally Exhausted","γδ T lymphocytes","NK T cells","CD56 Dim NK cells","CD56 Bright NK cells","Pre B cells","B cells","Pro B cells","Plasma Cells",  "Early Erythroids","Mid Erythroids","Late Erythroids","Blasts","Monocytes","Non Classical Monocytes","Pro Monocytes","cDC","pDC","HSPCs","Unidentified","Doublets")

ordered_timepoints= c("pre-transplant", "remission(3-6mo)","remission(>6mo)")

t_sum1$celltype = factor(t_sum1$celltype, levels = ordered_celltypes)
t_sum1$groups= factor(t_sum1$groups, levels = ordered_timepoints)

#Plot 
plot = ggplot(t_sum1, aes(x =  as.factor(groups), fill = celltype, group = celltype,
                          stratum = celltype, alluvium = celltype, 
                          y = freq, label = celltype)) +  geom_stratum() + geom_flow(stat = "alluvium") + geom_text(size=3, stat = "stratum", aes(label = after_stat(stratum))) + labs(x="timepoints")
plot


########## DGE ############
#Load the needed libraries
library(tidyverse)
library(ggplot2)
library(readxl)
library(data.table)
library(ggrepel)

#CD4 Memory cells in remission samples cohort  1 vs cohort 2
CD4M <- subset(T_Remission, subset =celltype %in% c("CD4 Memory"))


Idents(CD4M) <- CD4M$cohort
Analysis2<- FindMarkers(CD4M, ident.1 = "cohort1" , 
                        ident.2 = "cohort2",
                        logfc.threshold = 0.25,
                        test.use = "wilcox",
                        min.pct = 0.1)

data2 = Analysis2 %>% rownames_to_column("gene")

#Add a column to identify if it the gene is upregulated or down regulated according to adj p value and AvgLog2FoldChange
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
data2$diffexpressed[data2$avg_log2FC > 0.6 & data2$p_val_adj < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
data2$diffexpressed[data2$avg_log2FC < -0.6 & data2$p_val_adj < 0.05] <- "DOWN"
data2$gene[data2$diffexpressed != "NO"] <- data2$gene_symbol[data2$diffexpressed != "NO"]

#Add specific color pattern
mycolors <- c("blue", "red", "black")
names(mycolors) <- c("DOWN", "UP", "NO")

#Add this line to show p values those are 0 on the plot because when you take the log10 of it it makes it infinitive and makes it impossible to show them on the plot even changing the y axis limits
data2 = data2 %>% mutate(p_val_adj = ifelse(p_val_adj==0, 1e-300, p_val_adj))

p2 <- ggplot(data=data2, aes(x=avg_log2FC, y=-log10(p_val_adj), col=diffexpressed, label=gene)) +
  geom_point()+
  theme(aspect.ratio = 0.5,
        axis.text.x = element_text(angle = 45, vjust= 1, hjust = 1, size = 15, color = "black"),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 15, color = "black"),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        plot.title = element_text(size = 14, hjust = 0.5))+
  theme_bw()+
  geom_text_repel()+
  geom_vline(xintercept=c(-0.6, 0.6), col="black") +
  geom_hline(yintercept=-log10(0.05), col="black")
p2 

#C8 cells in remission samples cohort  1 vs cohort 2

CD8_remission <- subset(T_Remission, subset =celltype %in% c("CD8 Effector Memory", "CD8 Central Memory", "CD8 Naïve", "CD8 Terminally Exhausted"," γδ T lymphocyte"))


Idents(CD8_remission) <- CD8_remission$cohort
Analysis2<- FindMarkers(CD8_remission, ident.1 = "cohort1" , 
                        ident.2 = "cohort2",
                        logfc.threshold = 0.25,
                        test.use = "wilcox",
                        min.pct = 0.1)

data2 = Analysis2 %>% rownames_to_column("gene")

#Add a column to identify if it the gene is upregulated or down regulated according to adj p value and AvgLog2FoldChange
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
data2$diffexpressed[data2$avg_log2FC > 0.6 & data2$p_val_adj < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
data2$diffexpressed[data2$avg_log2FC < -0.6 & data2$p_val_adj < 0.05] <- "DOWN"
data2$gene[data2$diffexpressed != "NO"] <- data2$gene_symbol[data2$diffexpressed != "NO"]

#Add specific color pattern
mycolors <- c("blue", "red", "black")
names(mycolors) <- c("DOWN", "UP", "NO")

#Add this line to show p values those are 0 on the plot because when you take the log10 of it it makes it infiniotive and makes it impossible to show them on the plot even changin the y axis limits
data2 = data2 %>% mutate(p_val_adj = ifelse(p_val_adj==0, 1e-300, p_val_adj))

p2 <- ggplot(data=data2, aes(x=avg_log2FC, y=-log10(p_val_adj), col=diffexpressed, label=gene)) +
  geom_point()+
  theme(aspect.ratio = 0.5,
        axis.text.x = element_text(angle = 45, vjust= 1, hjust = 1, size = 15, color = "black"),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 15, color = "black"),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        plot.title = element_text(size = 14, hjust = 0.5))+
  theme_bw()+
  geom_text_repel()+
  geom_vline(xintercept=c(-0.6, 0.6), col="black") +
  geom_hline(yintercept=-log10(0.05), col="black")
p2


#C8 cells in remission samples remission vs relapse samples
CD8 <- subset(tc_diet, subset =celltype %in% c("CD8 Effector Memory", "CD8 Central Memory", "CD8 Naïve", "CD8 Terminally Exhausted","у8 T lymphocytes"))

Idents(CD8) <- CD8$status
Analysis1<- FindMarkers(CD8, ident.1 = "remission" , 
                        ident.2 = "relapse",
                        logfc.threshold = 0.25,
                        test.use = "wilcox",
                        min.pct = 0.1)

data1 = Analysis1 %>% rownames_to_column("gene")
#if you exporting the dataset use this line of code if not skip it 
write.csv(Analysis1_tibble, file = "~/Analysis1.csv", quote = F, append = F, row.names = F)

#Add a column to identify if it the gene is upregulated or down regulated according to adj p value and AvgLog2FoldChange
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
data1$diffexpressed[data1$avg_log2FC > 0.6 & data1$p_val_adj < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
data1$diffexpressed[data1$avg_log2FC < -0.6 & data1$p_val_adj < 0.05] <- "DOWN"
data1$gene[data1$diffexpressed != "NO"] <- data1$gene_symbol[data1$diffexpressed != "NO"]

#Add specific color pattern
mycolors <- c("blue", "red", "black")
names(mycolors) <- c("DOWN", "UP", "NO")

#Add this line to show p values those are 0 on the plot because when you take the log10 of it it makes it infiniotive and makes it impossible to show them on the plot even changin the y axis limits
data1 = data1 %>% mutate(p_val_adj = ifelse(p_val_adj==0, 1e-300, p_val_adj))
p1 <- ggplot() + 
  geom_point(data=data1, aes(x=avg_log2FC, y=-log10(p_val_adj), col=diffexpressed)) + 
  theme(aspect.ratio = 0.5,
        axis.text.x = element_text(angle = 45, vjust= 1, hjust = 1, size = 15, color = "black"),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 15, color = "black"),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        plot.title = element_text(size = 14, hjust = 0.5))+
  theme_bw()+
  #upcoming line we filtered the genes since we had alot next time you can adjust the number and made overlaps infinite
  geom_text_repel(data=data1 %>% filter(-log10(p_val_adj) > 100), aes(x=avg_log2FC, y=-log10(p_val_adj), label=gene), max.overlaps = Inf)+
  geom_vline(xintercept=c(-0.6, 0.6), col="black") +
  geom_hline(yintercept=-log10(0.05), col="black")
p1

p2 <- ggplot(data=data1, aes(x=avg_log2FC, y=-log10(p_val_adj), col=diffexpressed, label=gene)) +
  geom_point()+
theme(aspect.ratio = 0.5,
      axis.text.x = element_text(angle = 45, vjust= 1, hjust = 1, size = 15, color = "black"),
      axis.title.x = element_blank(),
      axis.text.y = element_text(size = 12),
      axis.title.y = element_text(size = 15, color = "black"),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 12),
      plot.title = element_text(size = 14, hjust = 0.5))+
  theme_bw()+
  geom_text_repel()+
  geom_vline(xintercept=c(-0.6, 0.6), col="black") +
  geom_hline(yintercept=-log10(0.05), col="black")
p2

########Subset the blasts populations#######

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

######## Pseudobulk DGE analysis #######
#Load the needed libraries
library(tidyverse)
library(Seurat)
library(ggplot2)
library(harmony)
library(randomcoloR)
library(RColorBrewer)
library(limma)
library(readxl)
library(cowplot)
library(SingleCellExperiment)
library(ggpubr)
library(data.table)
library(pheatmap)
library(apeglm)
library(png)
library(DESeq2)
library(Matrix.utils)
library(ggforce)
library(cloudml)

View(ex@meta.data)
sce = as.SingleCellExperiment(ex)
celltype_names = unique(sce$celltype)
length(celltype_names)

sce$sample_id = paste0(sce$id)
head(colData(sce))

# Subset metadata to include only the variables you want to aggregate across (here, we want to aggregate by sample and by celltype)

groups = colData(sce)[, c("celltype", "sample_id")]
head(groups)

# Aggregate across cluster-sample groups
# transposing row/columns to have cell_ids as row names matching those of groups
aggr_counts = aggregate.Matrix(t(counts(sce)), 
                               groupings = groups, fun = "sum") 
# Explore output matrix
class(aggr_counts)
dim(aggr_counts)
aggr_counts[1:6, 1:6]

# Transpose aggregated matrix to have genes as rows and samples as columns
aggr_counts = t(aggr_counts)
aggr_counts[1:6, 1:6]


tstrsplit(colnames(aggr_counts), "_") %>% str()
head(colnames(aggr_counts), n = 10)
head(tstrsplit(colnames(aggr_counts), "_")[[1]], n = 10)

#test fot b cells
#b_cell_idx <- which(tstrsplit(colnames(aggr_counts), "_")[[1]] == "B cells")
#b_cell_idx
#colnames(aggr_counts)[b_cell_idx]
#aggr_counts[1:10, b_cell_idx]

celltype_names
## Initiate empty list
counts_ls = list()

#put in a loop
for (i in 1:length(celltype_names)) {
  
  ## Extract indexes of columns in the global matrix that match a given cluster
  column_idx = which(tstrsplit(colnames(aggr_counts), "_")[[1]] == celltype_names[i])
  
  ## Store corresponding sub-matrix as one element of a list
  counts_ls[[i]] = aggr_counts[, column_idx]
  names(counts_ls)[i] = celltype_names[i]
  
}

# Explore the different components of the list
str(counts_ls)

head(colData(sce))

# Extract sample-level variables
metadata = colData(sce) %>% 
  as.data.frame() %>% 
  dplyr::select(patient_identity, status, cohort, sample_id)

metadata = metadata[!duplicated(metadata), ]

dim(metadata)
head(metadata)

rownames(metadata) = metadata$sample_id
head(metadata)

# Number of cells per sample and cluster
t <- table(colData(sce)$sample_id,
           colData(sce)$celltype)
t[1:6, 1:6]

# Creating metadata list
colnames(counts_ls)
## Initiate empty list
metadata_ls <- list()

for (i in 1:length(counts_ls)) {
  
  ## Initiate a data frame for cluster i with one row per sample (matching column names in the counts matrix)
  df <- data.frame(celltype_sample_id = colnames(counts_ls[[i]]))
  
  ## Use tstrsplit() to separate cluster (cell type) and sample IDs
  df$celltype <- tstrsplit(df$celltype_sample_id, "_")[[1]]
  df$sample_id  <- tstrsplit(df$celltype_sample_id, "_")[[2]]
  
  
  ## Retrieve cell count information for this cluster from global cell count table
  idx <- which(colnames(t) == unique(df$celltype))
  cell_counts <- t[, idx]
  
  ## Remove samples with zero cell contributing to the cluster
  cell_counts <- cell_counts[cell_counts > 0]
  
  ## Match order of cell_counts and sample_ids
  sample_order <- match(df$sample_id, names(cell_counts))
  cell_counts <- cell_counts[sample_order]
  
  ## Append cell_counts to data frame
  df$cell_count <- cell_counts
  
  ## Join data frame (capturing metadata specific to cluster) to generic metadata
  df <- plyr::join(df, metadata, 
                   by = intersect(names(df), names(metadata)))
  
  ## Update rownames of metadata to match colnames of count matrix, as needed later for DE
  rownames(df) <- df$cluster_sample_id
  
  ## Store complete metadata for cluster i in list
  metadata_ls[[i]] <- df
  names(metadata_ls)[i] <- unique(df$celltype)
  
}

# Explore the different components of the list
str(metadata_ls)


####### DEG #####
# Select cell type of interest
celltype_names

# Double-check that both lists have same names
all(names(counts_ls) == names(metadata_ls))

##
#From the website, run this when you are running one celltype ----- test run#
idx <- which(names(counts_ls) == "CD8 Terminally Exhausted")
cluster_counts <- counts_ls[[idx]]
cluster_metadata <- metadata_ls[[idx]]

# Check contents of extracted objects
cluster_counts[1:6, 1:6]
head(cluster_metadata)

#run this code to correct rownames when running for only one cell type
cluster_metadata = cluster_metadata %>% `row.names<-`(cluster_metadata$celltype_sample_id)

# Check matching of matrix columns and metadata rows
all(colnames(cluster_counts) == rownames(cluster_metadata))

#----------#
#Ksenia's code- for the loop of all the celltypes
total_res_table_thres = data.frame()

# Iterate over all celltypes
for (celltype in unique(sce$celltype)) {
  print(celltype)
  
  # Only run DE if cell counts for both conditions exceed 99
  cell_counts = metadata_ls[[celltype]] %>% as_tibble() %>% group_by(cohort) %>% summarize(n=sum(cell_count)) %>% mutate(ind = n >=100) %>% pull(ind) %>% sum()
  if (cell_counts < 2) {
    print(paste0(celltype, ' doesn\'t have enough cells, skip'))
    next
  }
  idx = which(names(counts_ls) == celltype)
  cluster_counts = counts_ls[[idx]]
  cluster_metadata = metadata_ls[[idx]]
  
  
# Create DESeq2 object 
dds = DESeqDataSetFromMatrix(cluster_counts, 
                             colData = cluster_metadata, 
                             design =  ~ cohort)

# Transform counts for data visualization
rld = rlog(dds, blind=TRUE)

# Plot PCA
p1 = DESeq2::plotPCA(rld, ntop = 500, intgroup = "cohort") + ggrepel::geom_text_repel(aes(label = name)) + theme_pubr()

DESeq2::plotPCA(rld, ntop = 500, intgroup = "cell_count")

# Extract the rlog matrix from the object and compute pairwise correlation values
rld_mat = assay(rld)
rld_cor = cor(rld_mat)

annotation = cluster_metadata[, c("cohort"), drop=F] %>% rownames_to_column("temp") %>% mutate(temp = sub("-.\\.", ".", temp)) %>% unique() 
rownames(annotation) = annotation$temp

#Plot heatmap
#p2 = pheatmap(rld_cor, annotation = annotation %>% dplyr::select(-temp))

 # Run DESeq2 differential expression analysis
 dds <- DESeq(dds, quiet = T)

 #Plot dispersion estimates
 #plotDispEsts(dds)
 
 # Check the coefficients for the comparison
 #resultsNames(dds)

 # Generate results object
 res <- results(dds, 
                name = "cohort_cohort2_vs_cohort1",
                alpha = 0.05)
 
 # Shrink the log2 fold changes to be more appropriate using the apeglm method - should cite [paper]() when using this method
 res = lfcShrink(dds, 
                 coef = "cohort_cohort2_vs_cohort1",
                 res=res,
                 type = "ashr",
                 quiet = T)

# Shrink the log2 fold changes to be more appropriate using the apeglm method - should cite [paper]() when using this method
#res <- lfcShrink(dds, 
#                coef = "cohort_cohort2_vs_cohort1",
#                res=res,
#               type = "apeglm") 
 
 # generate the results table for all of our genes, ordered by adjusted p-value 
 # Turn the DESeq2 results object into a tibble for use with tidyverse functions
 res_tbl = res %>%
   data.frame() %>%
   rownames_to_column(var = "gene") %>%
   as_tibble() %>%
   arrange(padj)

# Check results output
#res_tbl

# Write all results to file
#write.csv(res_tbl,
#         paste0("~/", unique(cluster_metadata$cluster_id), "_", 
#                levels(cluster_metadata$group_id)[2], "_vs_", levels(cluster_metadata$group_id)[1], "_all_genes.csv"),
#         quote = FALSE, 
#         row.names = FALSE)

# Set thresholds
padj_cutoff <- 0.05
# Subset the significant results
sig_res <- dplyr::filter(res_tbl, padj < padj_cutoff) %>%
  dplyr::arrange(padj)

#sig_res
#Ksenia: My code crashes when there is less than two DE genes and I am not motivated to fix that yet
if (nrow(sig_res) < 2) {
  print(paste0(celltype, ' doesn\'t have enough DEGs, skip'))
  next
}

## Extract normalized counts from dds object
normalized_counts = counts(dds, normalized = TRUE)

## Extract top 20 DEG from resLFC (make sure to order by padj)
top20_sig_genes <- sig_res %>%
  dplyr::arrange(padj) %>%
  dplyr::pull(gene) %>%
  head(n = 20)

## Extract matching normalized count values from matrix
top20_sig_counts <- normalized_counts[rownames(normalized_counts) %in% top20_sig_genes, ]
#top20_sig_counts

## Convert wide matrix to long data frame for ggplot2
top20_sig_df <- data.frame(top20_sig_counts)
top20_sig_df$gene <- rownames(top20_sig_counts)

#This part crashes if you don't have more than one single gene that is log fold different 
top20_sig_df <- melt(setDT(top20_sig_df), 
                     id.vars = c("gene"),
                     variable.name = "celltype_sample_id") %>% 
data.frame()

# Set thresholds
log2fc_cutoff <- 0.58

# Count significantly up/down genes above threshold
n_sig_up <- dplyr::filter(sig_res, log2FoldChange >= log2fc_cutoff) %>% 
nrow()
n_sig_dn <- dplyr::filter(sig_res, log2FoldChange <= -log2fc_cutoff) %>% 
nrow()

## Replace "." by " " in cluster_sample_id variable (melt() introduced the ".")
top20_sig_df$celltype_sample_id <- gsub("\\.", " ", top20_sig_df$celltype_sample_id)
top20_sig_df

## Join counts data frame with metadata
top20_sig_df <- plyr::join(top20_sig_df, as.data.frame(colData(dds)),
                           by = "celltype_sample_id")
top20_sig_df

## Generate plot
p3= ggplot(top20_sig_df, aes(y = value, x = cohort, col = cohort)) +
  geom_jitter(height = 0, width = 0.15) +
  scale_y_continuous(trans = 'log10') +
  ylab("log10 of normalized expression level") +
  xlab("condition") +
  ggtitle("Top 20 Significant DE Genes") +
  theme(plot.title = element_text(hjust = 0.5)) +
  facet_wrap(~ gene)
p3

# Heatmap

## Extract normalized counts for significant genes only
sig_counts <- normalized_counts[rownames(normalized_counts) %in% sig_res$gene, ]

sig_counts
## Set a color-blind friendly palette
heat_colors <- rev(brewer.pal(11, "PuOr"))

## Run pheatmap using the metadata data frame for the annotation
p4 = pheatmap(sig_counts, 
         color = heat_colors, 
         cluster_rows = TRUE, 
         show_rownames = FALSE,
         annotation = cluster_metadata[, c("sample_id", "celltype")], 
         border_color = NA, 
         fontsize = 10, 
         scale = "row", 
         fontsize_row = 10, 
         height = 40)
         #filename = paste0("", sub('/', '_', celltype), ".pdf" ))  

p4
# Volcano plot
log2fc_cutoff = 0.26 # (20% increase)
res_table_thres = res_tbl[!is.na(res_tbl$padj), ] %>% 
  mutate(threshold = padj < padj_cutoff & abs(log2FoldChange) >= log2fc_cutoff) 

top_genes = top20_sig_df$gene %>% unique()
top_genes
## Generate plot
p5 = ggplot() +
  geom_point(data=res_table_thres, aes(x = log2FoldChange, y = -log10(padj), colour = threshold)) +
  ggrepel::geom_text_repel(data=res_table_thres %>% filter(gene %in% top_genes), aes(x = log2FoldChange, y = -log10(padj), label = gene)) +
  xlab("log2 fold change") +
  ylab("-log10 adjusted p-value") +
  scale_color_manual(values = c("grey60", "red3")) +
  theme_pubr() + 
  ggtitle("Cohort1 vs Cohort2 DE")
p5

write.csv(x = res_table_thres, file = "~/CD8  Memory_tib.csv")



res_table_thres$celltype = celltype
total_res_table_thres = rbind(total_res_table_thres, res_table_thres)
p = ggarrange(ggarrange(p1,p5, ncol=1, widths = c(0.3,0.7), align = "none"),p3, align = "none")
save_plot(paste0("", sub('/', '_', celltype), ".summary.pdf"), p, base_height = 20, base_width = 20)

}

### Write all results to file
name = "all_celltypes"
write.csv(x = total_res_table_thres, file = paste0(name, ".pseudobulk_DE_res.csv"), quote = F, row.names = F)


######## GSEA ######
seu_diet_merged@active.ident <- seu_diet_merged$cohort
Idents(seu_diet_merged) <- seu_diet_merged$cohort
Analysis1<- FindMarkers(seu_diet_merged, ident.1 = "cohort1" , 
                        ident.2 = "cohort2",
                        logfc.threshold = 0.25,
                        test.use = "wilcox",
                        min.pct = 0.1)

Analysis1 <- Analysis1[order(Analysis1$avg_log2FC,decreasing =T),]
Genes <- rownames(Analysis1)

gene_rank <- Analysis1$avg_log2FC
names(gene_rank) <- Genes

h_gene_sets = msigdbr(species = "human", category = "H")
h_gene_sets = h_gene_sets %>% split(x = .$gene_symbol, f = .$gs_name)

gsea_out <-fgsea(h_gene_sets, stats = gene_rank)

# Get the data in the correct format

gsea_data <- data.frame(matrix(NA, ncol = 1, nrow = 42))
rownames(gsea_data) <- gsea_out$pathway
for(row in 1:nrow(gsea_data)){
  pvalue <- gsea_out$padj[row]
  if(pvalue < 0.05){
    gsea_data[row,] <- gsea_out$NES[row]
  }
  else(gsea_data[row,] <- 0)
}


# Make heatmap
# Set up for the heatmap
myColor=colorRampPalette(c("blue", "white", "red"))(100)
vals= gsea_data$matrix.NA..ncol...1..nrow...42.
myBreaks <- c(seq(min(vals), 0, length.out=ceiling(100/2) + 1), 0,
              seq(max(vals)/100, max(vals), length.out=floor(100/2)))
myBreaks2 = myBreaks + seq_along(myBreaks) * .Machine$double.eps
# Make the heatmap
install.packages("pheatmap")
library(pheatmap)
pheatmap(mat =gsea_data,
         color=myColor,
         cluster_cols = F,
         show_colnames = T,
         fontsize_row = 5,
         width=15,
         cellwidth = 10,
         cellheight = 10,
         breaks=myBreaks2,
         height = 15)


####### soupor cell ######

#calculate variant covarage 
#Load the matrices for 2737

M1 <- readMM("~/ref.mtx")
M2 <- readMM("~/alt.mtx")

#Load the matrices for 1013
M1 <- readMM("~/souporcell_1013_ref.mtx")
M2 <- readMM("~/souporcell_1013_alt.mtx")

df1 <- as.data.frame(summary(M1))
df2 <- as.data.frame(summary(M2))

#add a new column in to df1 sum of reference and alternative matrices
df1$cov <- df1$x + df2$x
#select only covarege 
df1 <- subset(df1, select = c("i","cov"))
#put into a table
table(df1$cov)
#Visualize it
p1 <- ggplot(df1 %>% filter( cov<30 ), aes(x=cov)) +
  geom_histogram( binwidth=3, fill="#69b3a2", color="#e9ecef", alpha=0.9)

p1

###adding souporcell annotations###

#to add a color palette
palette <- distinctColorPalette(k = 34)
pie(rep(1, 34), col = palette) 

#view the metadata
View(seu_diet@meta.data)

#load souporcell output tsv file that contains the assignments from the souporcell
df1 <- read_tsv("~/souporcell_Pt2_clusters.tsv")

#these lines only apply if you are running single bam file and we are trying yo wrangle the thing similar to seu mata data but you can also remove this from seu object as well
df1$barcode <- paste0(df1$barcode,"_27")

df2 <- seu_diet@meta.data
df2 = df2 %>% rownames_to_column("barcode")

assignment = df1$assignment
names(assignment) <- df1$barcode

seu_diet<- AddMetaData(object = seu_diet, metadata = assignment, col.name = "souporcell")

T2737 <- subset(x = seu_diet, subset = orig.ident %in% c("2737_MNC"))

##### You need to export barcodes to run souporcell.
df2 <- seu_diet_merged@meta.data
df2 = df2 %>% rownames_to_column("barcode")
P1 <- subset(x = df2, subset = orig.ident %in% c("2446_MNC", "25802_MNC", "25802_CD3", "2645_MNC"))
P2 <- subset(x = df2, subset = orig.ident %in% c("1972_MNC", "2220_MNC", "2621_MNC", "2621_CD3","1972_CD3"))

P3 <- subset(x = df2, subset = orig.ident %in% c("9185_MNC","9185_CD3"))
P4 <- subset(x = df2, subset = orig.ident %in% c("2599_MNC","2599_CD3"))
P7 <- subset(x = df2, subset = orig.ident %in% c("2518_MNC","2518_CD3"))

P5 <- subset(x = df2, subset = orig.ident %in% c("9596_MNC", "25809_MNC", "2737_MNC", "9596_CD3","2737_CD3"))
P6 <- subset(x = df2, subset = orig.ident %in% c("2379_MNC", "2434_MNC","2379_CD3"))
P8 <- subset(x = df2, subset = orig.ident %in% c("4618_MNC", "6174_MNC","9931_MNC","1953_MNC","6174_CD3","9931_CD3","1953_CD3"))


P9 <-  subset(x = df2, subset = orig.ident %in% c("1677_MNC", "1732_MNC","1811_MNC","1677_CD3","1811_CD3"))
P10 <- subset(x = df2, subset = orig.ident %in% c("1195_MNC", "1285_MNC","1347_MNC","1347_CD3"))
P11 <- subset(x = df2, subset = orig.ident %in% c("5641_MNC", "6244_MNC","6244_CD3"))
P12 <-  subset(x = df2, subset = orig.ident %in% c("9355_MNC", "1013_MNC")) 

barcode4 <- P4$barcode
#remove the sample identifier after 
barcode4 = gsub("_.*", "", barcode4)
#check the barcodes 
length(barcode4)
#exclude the doublets that across the different samples
barcode4_unique = as.data.frame(table(barcode4)) %>% filter(Freq == 1) %>% pull(barcode4)
print(barcode4_unique, row.names = F)
#check again to see how many of them dropped
length(barcode4_unique)

(length(barcode4)-length(barcode4_unique))/length(barcode4)

#export as tsv file since souporcell requires that
write_tsv((as.data.frame(barcode4_unique)), file = "~/pt4.barcodes.tsv", col_names = FALSE)

##these lines apply when using merged bam file
P1 <- subset(x = seu_diet_merged, subset = orig.ident %in% c("2446_MNC","25802_MNC","2645_MNC","25802_CD3"))
P2 <- subset(x = seu_diet_merged, subset = orig.ident %in% c("1972_MNC","2220_MNC","2621_MNC","1972_CD3","2621_CD3"))


P3 <- subset(x = seu_diet_merged, subset = orig.ident %in% c("9185_MNC","9185_CD3"))
P4 <- subset(x = seu_diet_merged, subset = orig.ident %in% c("2599_MNC","2599_CD3"))
P7 <- subset(x = seu_diet_merged, subset = orig.ident %in% c("2518_MNC","2518_CD3"))


P5 <- subset(x = seu_diet_merged, subset = orig.ident %in% c("9596_MNC","25809_MNC","2737_MNC","9596_CD3","2737_CD3"))
P6 <- subset(x = seu_diet_merged, subset = orig.ident %in% c("2379_MNC","2434_MNC","2379_CD3"))
P8 <- subset(x = seu_diet_merged, subset = orig.ident %in% c("4618_MNC","6174_CD3","9931_CD3","1953_CD3","6174_MNC","9931_MNC","1953_MNC"))

P9 <- subset(x = seu_diet_merged, subset = orig.ident %in% c("1677_MNC","1677_CD3","1732_MNC","1811_MNC","1811_CD3"))
P10 <- subset(x = seu_diet_merged, subset = orig.ident %in% c("1195_MNC","1285_MNC","1347_MNC","1347_CD3"))
P11 <- subset(x = seu_diet_merged, subset = orig.ident %in% c("5641_MNC","6244_MNC","6244_CD3"))
P12 <- subset(x = seu_diet_merged, subset = orig.ident %in% c("9355_MNC","1013_MNC"))


#load souporcell output tsv file that contains the assignments from the souporcell
Df1 <- read_tsv("~/souporcell_Pt4_clusters.tsv")

#wrangle the df
P4 = as.data.frame(P4@meta.data) %>% rownames_to_column("barcode")
P4$barcode <- gsub("_.*", "", P4$barcode)


#merge 2 data frame  by using the joint column, it will add to the first df you write
df4 <- P4%>% left_join(Df1, by="barcode")

library(janitor)

#Check the actual cell numbers with total cell numbers at the end
df4 %>%
     tabyl(celltype, assignment,id , show_missing_levels = FALSE) %>%
     adorn_totals("row")

#remove the doublet assignments
df4 <- subset(df4, subset = assignment %in% c("0","1")) 

#rename the assignment for better read
df4$assignment <- gsub("0", "donor", df4$assignment)
df4$assignment <- gsub("1", "host", df4$assignment) 

df4$status <- df4$status.x
df4 <- df4 %>% 
  dplyr::select(orig.ident, celltype, library_type, patient_identity, cell, id, status, assignment, cohort)

df4$cohort <- gsub("cohort1", "Remission cohort", df4$cohort)


combined_df1 <- bind_rows(df1,df2)
write_csv(combined_df1, "~/cohort1_souporcell.csv")
combined_df1 <- read_csv(file = "~/cohort1_souporcell.csv")

 
combined_df2 <- bind_rows(df5,df6,df8)
write_csv(combined_df2, "~/cohort2_souporcell.csv")
combined_df2 <- read_csv(file = "~/cohort2_souporcell.csv")

#After long discussion with Peter we have decided to assign all of the Patient 3 and Patient 7 samples as donor since the chimerism data was higher than 95 percent and we did not have a pre-transplant samples as we did in others which caused souporcell to give inaccurate results.

P7$assignment <- c("donor")
P3$assignment <- c("donor")


df1 <- 
combined_df1 %>%
dplyr::select(orig.ident, celltype, library_type, patient_identity, cell, id, status.x, assignment, cohort)

df1$status <- df1$status.x
df1 <- df1 %>% 
  dplyr::select(orig.ident, celltype, library_type, patient_identity, cell, id, status, assignment, cohort)
  
P3 = as.data.frame(P3@meta.data) %>% rownames_to_column("barcode")
View(P3)

df3 <- P3 %>% 
  dplyr::select(orig.ident, celltype, library_type, patient_identity, cell, id, status, assignment,cohort)

View(df3)

#add P3 to cohort 1 dataframe
combined_df1 <- bind_rows(df1,df3)

df2 <- 
  combined_df2 %>%
  dplyr::select(orig.ident, celltype, library_type, patient_identity, cell, id, status.x, assignment, cohort)

df2$status <- df2$status.x
df2 <- df2 %>% 
  dplyr::select(orig.ident, celltype, library_type, patient_identity, cell, id, status, assignment, cohort)

P7 = as.data.frame(P7@meta.data) %>% rownames_to_column("barcode")
View(P7)

df7 <- P7 %>% 
  dplyr::select(orig.ident, celltype, library_type, patient_identity, cell, id, status, assignment,cohort)

View(df7)

#add P7 to cohort 2 dataframe
combined_df2 <- bind_rows(df2,df7)


View(combined_df2)

combined_df2$cohort <- gsub("cohort2", "Relapse cohort", combined_df2$cohort) 
combined_df1$cohort <- gsub("cohort1", "Remission cohort", combined_df1$cohort)


combined_df3 <- bind_rows(df9,df10,df11,df12)
write_csv(combined_df3, "~/cohort3_souporcell.csv")
combined_df3 <- read_csv(file = "~/cohort3_souporcell.csv")



#load the combined dataframe
combined_df<- bind_rows(combined_df1,combined_df2)
write_csv(combined_df, "~/cohort1-2_souporcell.csv")
combined_df <- read_csv(file = "~/cohort1-2_souporcell.csv")

View(combined_df)

#Pool MNC and CD3 samples
combined_df <- combined_df%>%
  mutate(pt_timepoint = gsub("RemT", "Rem", id)) 

#select only T cells(no NK cells)
combined_df_T <- subset(combined_df, subset = celltype %in% c("CD4 Memory","CD8 Effector","CD8 Memory","CD4 Naïve","Treg","CD8 Naïve","CD8 Effector","CD4 Naïve","CD8 Terminally Exhausted","CD4 Memory","NK T cells","γδ T lymphocytes", "CD4 Naïve", "CD56 Bright NK cells", "CD56 Dim NK cells"))

#select only remission cells
combined_df_T_rem <- subset(combined_df_T, subset = status == "remission")
#select only 6mo remission cells and MNC libraries only to prevent any skewing
combined_df_T_rem6mo <- subset(combined_df_T_rem, subset = pt_timepoint %in% c("P01.1Rem","P01.2Rem","P02.1Rem","P04.1Rem","P05.1Rem","P06.1Rem","P07.1Rem","P08.1Rem"))

View(combined_df_T_rem6mo)

#select only donor cells
combined_df_T_rem6modonor <- subset(combined_df_T_rem6mo, subset = assignment == "donor")


#to export the table 
t1 <- combined_df_T_rem6modonor %>% 
  group_by(celltype, cohort, pt_timepoint, assignment) %>%
  summarize(n = n()) %>%
  group_by(pt_timepoint) %>%
  mutate(freq = n/sum(n)) 
  #add this line fot the average each patient sample
  #%>%
  #group_by(celltype, cohort, patient_identity) %>%
  #summarize(freq = mean(freq))



t2 <- combined_df  %>% 
  group_by(id,assignment) %>%
  summarize(n = n()) %>%
  mutate(freq = n/sum(n)) 

t3 <- combined_df3  %>% 
  group_by(id,assignment) %>%
  summarize(n = n()) %>%
  mutate(freq = n/sum(n)) 
head(t3)

write_csv(t2, "~/t2.csv")
write_csv(t3, "~/t3.csv")

head(t1)

#check if the each sample adds up the 1 
t1 %>% filter(pt_timepoint == "P07.1Rem") %>% pull(freq) %>% sum()


combined_df_T_rem6modonor %>% filter(grepl("NK", celltype)) %>% dim()
combined_df_T_rem6modonor %>% filter(grepl("CD56", celltype)) %>% dplyr::select(pt_timepoint) %>% table()

t1 %>%
  
  #filter(grepl("CD56", celltype)) %>%
  ggplot(aes(x=cohort, y=freq, color=cohort)) +
  geom_point(position = position_jitterdodge()) +
  facet_wrap(.~celltype) +
  stat_compare_means() +
  theme_pubr()
head(t1)


#export the table
 write_csv(t1, "~/t1.csv")
 

ordered_status = c("pre_transplant","remission","relapse")
meta_summary$status= factor(meta_summary$status, levels = ordered_status)

##Visualize##
#Alluvial graph
library(ggalluvial)

#change fill for geom_alluvium for to change the flow
ggplot(data = meta_summary,
       aes(axis1 = cohort , axis2 = assignment, axis3 = celltype,
           y = n_cells)) +
  scale_x_discrete(limits = c("cohort", "origin","celltype"), expand = c(.2, .05)) +
  geom_alluvium(aes(fill = celltype)) +
  geom_stratum(aes(fill=celltype), color="black") +
  geom_text(size= 3, stat = "stratum", aes(label =paste(after_stat(stratum)))) +
  scale_fill_manual(values = palette) + 
  ggtitle("Post-transplant 3-6 months remission samples")


## Plot showing cell proportions for  ##  
combined_donor6mo_T   %>% 
  select(id, celltype) %>% 
  group_by(id, celltype) %>% 
  summarize(n=n()) %>% 
  ggplot(aes(x= id, y=n, fill=celltype, label= n)) + 
  geom_bar(stat = "identity") + 
  geom_text(size = 3, position = position_stack(vjust = 0.5),colour = "white")+
  theme_pubr(base_size = 15) +
  theme(axis.text.x = element_text(angle=90, vjust = 0.5), legend.position = "right")

#actual cell numbers
combined_donor6mo_T %>% 
   select(id, assignment, celltype) %>% 
   group_by(celltype, assignment) %>% 
   summarize(n=n()) %>% 
   ggplot(aes(x= celltype, y=n, fill=assignment, label= n)) + 
   geom_bar(stat = "identity") + 
   geom_text(size = 3, position = position_stack(vjust = 0.5),colour = "white")+
   theme_pubr(base_size = 15) +
   scale_fill_manual(values = palette)+ 
  theme(axis.text.x = element_text(angle = 45, vjust= 1, hjust = 1, size = 10), legend.position = "right")

#to do t test 
library(ggforce)

 t1 %>%
  ggplot(aes(x=cohort, y=freq, color=cohort)) +
  geom_point(position = position_jitterdodge()) +
  facet_wrap(.~celltype) +
  stat_compare_means()

 
  t1 %>%
   filter(grepl("CD56", celltype)) %>%
   ggplot(aes(x=cohort, y=freq, color=cohort)) +
   geom_point(position = position_jitterdodge()) +
   facet_wrap(.~celltype) +
   stat_compare_means() +
    theme_pubr()
