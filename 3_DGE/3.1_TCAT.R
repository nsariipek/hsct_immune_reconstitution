# Peter van Galen, 250408
# Run TCAT and evaluate the algorithm's cell type labels

# Load libraries
library(Seurat)
library(tidyverse)
library(scattermore)
library(Matrix)
library(R.utils)

# Set working directory
#VM: setwd("~/TP53_ImmuneEscape/3_DGE/")
#Local: setwd("~/DropboxMGB/Projects/ImmuneEscapeTP53/TP53_ImmuneEscape/3_DGE")

# Delete environment variables & load favorite function
rm(list=ls())
cutf <- function(x, f=1, d="/") sapply(strsplit(x, d), function(i) paste(i[f], collapse=d))

# Load T cell data. This is available on Dropbox but not GitHub
seu_T <- readRDS("../AuxiliaryFiles/250128_Tcell_subset.rds")

# Load colors from 2.3_PvG-Colors.R
celltype_colors_df <- read.table("../celltype_colors.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE, comment.char = "")
celltype_colors <- setNames(celltype_colors_df$color, celltype_colors_df$celltype)

# Check data visually
seu_T@meta.data %>% sample_frac(1) %>%
          ggplot(aes(x = UMAP_TNK_1, y = UMAP_TNK_2, color = celltype)) +
          geom_scattermore(pointsize = 16, pixels = c(4096, 4096)) +
          scale_color_manual(values = celltype_colors) +
          theme_bw() +
          theme(aspect.ratio = 1,
                panel.grid = element_blank()) +
          guides(color = guide_legend(override.aes = list(size = 3)))

# Next, save Seurat count matrix to run TCAT (see https://github.com/immunogenomics/starCAT/blob/main/Examples/starCAT_vignette_R.ipynb)

# Extract T cell counts
counts <- LayerData(seu_T, assay = "RNA", layer = "counts")
counts[80:90, 1:10]

# Save counts matrix
dir.create("AuxiliaryFiles")
writeMM(counts, file = "AuxiliaryFiles/matrix.mtx")
gzip("AuxiliaryFiles/matrix.mtx", overwrite = T, remove = T)

# Save cell barcodes
barcodes <- colnames(counts)
write_delim(as.data.frame(barcodes),  file = "AuxiliaryFiles/barcodes.tsv", col_names = FALSE)
gzip("AuxiliaryFiles/barcodes.tsv", overwrite = T, remove = T)

# Save feature names
gene_names <- rownames(counts)
features <- data.frame("gene_id" = gene_names, "gene_name" = gene_names,type = "Gene Expression")
write_delim(as.data.frame(features),delim = "\t", file = "AuxiliaryFiles/features.tsv", col_names = FALSE)
gzip("AuxiliaryFiles/features.tsv", overwrite = T, remove = T)

# Run 3.2_PvG-TCAT.sh to generate results (this was done on a Google Cloud Platform virtual machine, bash). This will create a folder 3.1_starCAT with the results
#cd /home/unix/vangalen/TP53_ImmuneEscape/3_DGE
#mkdir 3.1_starCAT
#./3.1_starCAT.sh
#rm -r cache

# Zip the results to enable syncing over GitHub
gzip("3.1_starCAT/starCAT.rf_usage_normalized.txt", overwrite = T, remove = T)
gzip("3.1_starCAT/starCAT.scores.txt", overwrite = T, remove = T)

# The rest of the script compares TCAT cell annotations to our own and is not used for the paper. Start with loading scores.
scores_tib <- read_tsv("3.1_starCAT/starCAT.scores.txt.gz") %>% rename("cell" = "...1")

# Compare cell type annotations with Multinomial_Label from scores
metadata_tib <- as_tibble(seu_T@meta.data, rownames = "cell")
metadata_tib <- left_join(metadata_tib, scores_tib)

# Wrangle and assign colors
metadata_tib$Multinomial_Label <- factor(metadata_tib$Multinomial_Label, levels = c("CD4_Naive", "CD4_EM", "CD4_CM", "Treg",
      "CD8_Naive", "CD8_EM", "CD8_CM", "CD8_TEMRA", "gdT", "MAIT"))
TCAT_Label_colors <- c("CD4_Naive" = "#466983FF", "CD4_EM" = "#D58F5CFF", "CD4_CM" = "#C75127FF", "Treg" = "#FFC20AFF",
      "CD8_Naive" = "#33CC00FF", "CD8_EM" = "#612A79FF", "CD8_CM" = "#0099CCFF", "CD8_TEMRA" = "#CE3D32FF",
      "gdT" = "#D595A7FF", "MAIT" = "#0A47FFFF")

# UMAP
metadata_tib %>% sample_frac(1) %>%
          ggplot(aes(x = UMAP_TNK_1, y = UMAP_TNK_2, color = Multinomial_Label)) +
          geom_scattermore(pointsize = 16, pixels = c(4096, 4096)) +
          scale_color_manual(values = TCAT_Label_colors) +
          theme_bw() +
          theme(aspect.ratio = 1,
                panel.grid = element_blank()) +
          guides(color = guide_legend(override.aes = list(size = 3)))

# Compare cell type numbers
p1 <- metadata_tib %>% count(celltype) %>%
      ggplot(aes(x = celltype, y = n, fill = celltype)) +
      geom_bar(stat = "identity") +
      scale_fill_manual(values = celltype_colors) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
p2 <- metadata_tib %>% count(Multinomial_Label) %>%
      ggplot(aes(x = Multinomial_Label, y = n, fill = Multinomial_Label)) +
      geom_bar(stat = "identity") +
      scale_fill_manual(values = TCAT_Label_colors) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
p1 + p2

# Heatmap
df <- data.frame(celltype = rep(1:9, each = 10), Multinomial_Label = rep(1:9, times = 10))
metadata_tib %>% select(celltype, Multinomial_Label) %>%
      group_by(celltype, Multinomial_Label) %>% count() %>%
      #group_by(celltype) %>% 
      #mutate(prop = n / sum(n)) %>%
      ggplot(aes(x = Multinomial_Label, y = celltype, fill = n)) +
      geom_tile(color = "grey", size = 0.25) +
      geom_tile(data = df[df$celltype == df$Multinomial_Label, ], color = "black", size = 0.5, fill = NA) +
      scale_fill_gradient(low = "white", high = "red") + 
      theme_minimal() +
      labs(x = "Multinomial Label", y = "Cell Type", fill = "Count") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), aspect.ratio = 9/10)

# Compare annotations side by side (stacked cells)
metadata_tib %>%
      select(cell, celltype, Multinomial_Label) %>%
      arrange(desc(celltype), desc(Multinomial_Label)) %>%
      pivot_longer(cols = -cell, names_to = "annotation_approach", values_to = "anno") %>%
      mutate(cell = factor(cell, levels = unique(cell))) %>%
      ggplot(aes(x = annotation_approach, y = cell, fill = anno)) +
      geom_tile() +
      scale_fill_manual(values = c(celltype_colors, TCAT_Label_colors)) +
      theme_minimal() +
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
