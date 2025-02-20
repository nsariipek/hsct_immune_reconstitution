# Peter van Galen, 250217
# Combine variant-level data with cell metadata to make paper figures

# Load libraries
library(tidyverse)
library(Seurat)
library(ComplexHeatmap)
library(scales)
library(magick)

# Set working directory (local). For Nurefsan:
setwd("/Users/dz855/Dropbox (Partners HealthCare)/ImmuneEscapeTP53/TP53_ImmuneEscape/5_Souporcell/")
# For Peter:
#setwd("~/DropboxMGB/Projects/ImmuneEscapeTP53/TP53_ImmuneEscape/5_Souporcell/")

# Delete environment variables & load favorite function
rm(list=ls())
cutf <- function(x, f=1, d="/") sapply(strsplit(x, d), function(i) paste(i[f], collapse=d))

# Current patient to analyze
pt <- "P01"

# Load variants from 5.2_Variant_heatmaps.R
variants_tib <- read_tsv(file = paste0("variants/", pt, "_variants.txt"))

# Load Seurat data & extract relevant metadata
seu <- readRDS("../AuxiliaryFiles/250128_seurat_annotated_final.rds")
wrangled_metadata <- as_tibble(subset(seu, patient_id == pt)@meta.data, rownames = "barcode") %>%
  mutate(barcode = cutf(barcode, d = "_", f = 3)) %>% select(barcode, sample_status, celltype)

# Load cell type colors from 2.3_PvG-Colors.R
celltype_colors_df <- read.table("../celltype_colors.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE, comment.char = "")
celltype_colors <- setNames(celltype_colors_df$color, celltype_colors_df$celltype)

# Join with variant level data
heatmap_info_tib <- select(variants_tib, barcode, var, genotype, assignment) %>%
  left_join(wrangled_metadata) %>%
  mutate(genotype = as.character(genotype),
         assignment = factor(assignment))

# Restructure data more similar to the desired heatmap
heatmap_df <- heatmap_info_tib %>% slice_sample(n = nrow(.)) %>%
  pivot_wider(id_cols = c(barcode, sample_status, assignment, celltype), names_from = var, values_from = genotype) %>%
  as.data.frame() %>%
  mutate(sample_status = factor(sample_status, levels = c("pre_transplant", "remission", "relapse"))) %>%
  arrange(sample_status, assignment, celltype)

# Heatmap annotation
celltype_colors <- celltype_colors[factor(unique(heatmap_df$celltype),
  levels = names(celltype_colors))]
col_anno <- columnAnnotation(sample_status = heatmap_df$sample_status,
                             souporcell_assignment = heatmap_df$assignment,
                             celltype = heatmap_df$celltype,
                             col = list(sample_status = c("pre_transplant" = "#4775FFFF", "remission" = "#99CC00", "relapse" = "#BA6338"),
                              souporcell_assignment = c("0" = "#3B1B53", "1" = "#F0E685"),
                              celltype = celltype_colors))

# Convert numeric matrix to character matrix with "0" and "1"
heatmap_mat <- heatmap_df %>%
  column_to_rownames("barcode") %>%
  select(-c("sample_status", "assignment", "celltype")) %>%
  t()

# Generate heatmap
h1 <- Heatmap(heatmap_mat,
              col = c("0" = "#56B4E9", "1" = "#E69F00"),
              na_col = "white",
              cluster_rows = F,
              cluster_columns = F,
              use_raster = T, raster_quality = 5,
              bottom_annotation = col_anno,
              row_names_gp = gpar(fontsize = 6),
              column_names_gp = gpar(fontsize = 0),
              column_title = paste(comma(ncol(heatmap_mat)), "cells"),
              row_title = "Variant",
              heatmap_legend_param = list(title = "genotype"),
              border = T)

# Check plot
#h1

# Save
pdf(paste0("5.3_", pt, "_IntegratedHeatmap.pdf"), width = 12, height = 5)
h1
dev.off()


