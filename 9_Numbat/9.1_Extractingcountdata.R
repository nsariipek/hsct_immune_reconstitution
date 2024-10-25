#Nurefsan Sariipek, 240921
 #We need UMI counts for the second part of Numbat. I'll extract the corresponding UMI counts, save them as an RDS file, and copy it to my Google bucket so that I can transfer it to the Broad cluster to run the second part of Numbat

# Load the needed libraries
library(tidyverse)
library(Seurat)
library(readr)
library(dplyr)
library(numbat)

# Load the saved seurat objects
seu_diet_merged <- readRDS("/Users/dz855/Dropbox (Partners HealthCare)/ImmuneEscapeTP53/RDS files/seu_diet_merged.rds")

# Select only 2737 since that's the only one we run the Numbat, not sure if this is neccessary
#t <- subset(seu_diet_merged, subset = orig.ident %in%c("2737_MNC","25809_MNC","9596_MNC"))

t <- subset(seu_diet_merged, subset = orig.ident %in% c("2737_MNC","25809_MNC","9596_MNC"))
t3 <- subset(seu_diet_merged, subset = orig.ident =="2737_MNC")
t2 <- subset(seu_diet_merged, subset = orig.ident =="25809_MNC")
t1 <- subset(seu_diet_merged, subset = orig.ident =="9596_MNC")

# To ensure consistency with the allele data, the sample identifier following -1 in the cell names must be removed, do this uing the lines down below
# Get the current cell barcodes from the Seurat object 
cell_barcodes <- colnames(t)
# Modify the barcodes to keep '-1' and remove everything after
new_cell_barcodes <- sub("-1.*", "-1", cell_barcodes)
# Ensure that cell barcodes are unique by adding a suffix where necessary
new_cell_barcodes <- make.unique(new_cell_barcodes)
# Assign the modified cell names directly to the Seurat object 't1'
colnames(t) <- new_cell_barcodes
# Verify if the cell names are updated correctly
head(colnames(t))

# Subtract the UMI counts
umi_counts <- t@assays$RNA@counts

# Save as RDS file(similar to the example)
saveRDS(umi_counts, "pt5_numbat_umi_counts.rds")

