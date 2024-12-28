# Nurefsan Sariipek, 240921
# We need UMI counts for the second part of Numbat. For each patient, I'll extract the UMI counts, save them as an RDS file, and copy it to my Google bucket so that I can transfer it to the Broad cluster to run the second part of Numbat.

# Load the needed libraries
library(tidyverse)
library(Seurat)
library(readr)
library(numbat)

# Empty environment
rm(list=ls())

# Set working directory
# For Nurefsan
setwd("/Users/dz855/Dropbox (Partners HealthCare)/ImmuneEscapeTP53/")
# For Peter
#setwd("~/DropboxMGB/Projects/ImmuneEscapeTP53/")

# Load the saved seurat objects
seu_diet_merged <- readRDS("RDS files/seu_diet_merged.rds")

# Load the saved dataframe that contains souporcell information for selecting only host counts
soc_combined_df <- read_csv("TP53_ImmuneEscape/5_Souporcell/results/cohort1-2_souporcell.csv")
# Same info for cohort 3 
soc_combined_df <-read_csv("TP53_ImmuneEscape/5_Souporcell/results/cohort3_souporcell.csv")

# Subset souporcell output. Alternatively, use orig.ident %in% c("2737_MNC","25809_MNC","9596_MNC")
soc_subset <- soc_combined_df %>% filter(orig.ident %in% c("1195_MNC","1285_MNC","1347_MNC"), assignment == "host")

# # This code down below is not neccessary
# seu_diet_merged$cell = paste0(gsub("_.*", "", rownames(seu_diet_merged@meta.data)), "_", seu_diet_merged$orig.ident)
# soc_subset$cell = paste0(soc_subset$cell, "_", soc_subset$orig.ident)
# # Subset Seurat object for desired cells
# seu_subset <- subset(seu_diet_merged, cell %in% soc_subset$cell)

seu_subset <- subset(seu_diet_merged, cells = soc_subset$cell)

# Subtract the UMI counts
umi_counts <- seu_subset@assays$RNA@counts

# Modify the barcodes to keep '-1' and remove everything after
colnames(umi_counts) <- sub("-1.*", "-1", colnames(umi_counts))

# Check length
dim(umi_counts)

######################
# # To create reference data - we did not end up using this
# # Subset only donor cells 
# soc_subset <- soc_combined_df %>% filter(assignment == "donor")
# # Subset Seurat object for desired cells
# seu_subset <- subset(seu_diet_merged, cells = soc_subset$cell)
# # count_mat is a gene x cell raw count matrices
# # cell_annot is a dataframe with columns "cell" and "group"
# count_mat= seu_subset@assays$RNA@counts
# cell_annot= soc_subset %>%
#   select(celltype, cell) %>% 
#   rename(group = celltype) %>%
#   mutate(group = factor(group))
# 
# # Create the reference 
# ref_internal = aggregate_counts(count_mat, cell_annot)
# # Save as RDS file (similar to the example)
# saveRDS(ref_internal, "NS_reference.rds")
##########################

# Save as RDS file (similar to the example)
saveRDS(umi_counts, "pt10_host_numbat_umi_counts.rds")
