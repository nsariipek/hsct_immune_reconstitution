# Peter van Galen, 250528
# Add various types of metadata to the annotated Seurat object and save to facilitate downstream analyses

library(tidyverse)
library(Seurat)

# Clear environment
rm(list=ls())

# Set working directory
setwd("~/DropboxMGB/Projects/ImmuneEscapeTP53/TP53_ImmuneEscape/2_Annotate-predict/")

# Load the saved seurat objects
seu <- readRDS("../AuxiliaryFiles/250426_Seurat_annotated.rds")

# Add TCAT annotations
usage_tib <- read_tsv("../3_DGE/3.1_starCAT/starCAT.scores.txt.gz") %>% rename("cell" = "...1")
tcat_celltypes_df <- data.frame(select(usage_tib, Multinomial_Label, Proliferation), row.names = usage_tib$cell)
colnames(tcat_celltypes_df) <- c("TCAT_Multinomial_Label", "TCAT_Proliferation")
seu <- AddMetaData(seu, tcat_celltypes_df)
seu$TCAT_Multinomial_Label <- factor(seu$TCAT_Multinomial_Label, levels = c("CD4_Naive", "CD4_CM", "CD4_EM", "Treg", "CD8_Naive", "CD8_CM", "CD8_EM",  "CD8_TEMRA", "MAIT", "gdT"))

# Add TCR calls
tcr_calls <- read_csv("../7_TCR_Diversity/7.1_TCR_calls.csv.gz")
tcr_calls <- column_to_rownames(tcr_calls, var = "cell")
seu <- AddMetaData(seu, tcr_calls)

# Add souporcell results
souporcell_assignments <- read_csv("../6_Souporcell/6.2_Souporcell_assignments.csv.gz")
df_meta <- souporcell_assignments %>% select(cell, origin) %>% column_to_rownames("cell")
colnames(df_meta) <- "souporcell_origin"
seu <- AddMetaData(seu, metadata = df_meta)
seu$souporcell_origin <- factor(seu$souporcell_origin, levels = c("recipient", "donor"))

# Add Numbat data
numbat_calls <- read_csv("../8_Numbat/8.4_Numbat_calls.csv")
numbat_calls <- column_to_rownames(numbat_calls, var = "cell")
colnames(numbat_calls) <- "numbat_compartment"
seu <- AddMetaData(seu, numbat_calls)
seu$numbat_compartment <- factor(seu$numbat_compartment, levels = c("normal", "tumor"))

# View & save final object
seu@meta.data %>% head
saveRDS(seu, "../AuxiliaryFiles/250528_Seurat_complete.rds")


