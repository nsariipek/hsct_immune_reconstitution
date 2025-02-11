# Peter van Galen, 250206
# Run TCAT and evaluate the algorithm's cell type labels

library(Seurat)
library(tidyverse)
library(scattermore)
library(Matrix)
library(R.utils)

# Save updated Seurat object to home directory in bash:
#cd ~
#gsutil cp gs://fc-3783b423-62ac-4c69-8c2f-98cb0ee4503b/250128_seurat_annotated_final.rds .
#gsutil cp gs://fc-3783b423-62ac-4c69-8c2f-98cb0ee4503b/250128_Tcell_subset.rds .
#gsutil cp gs://fc-3783b423-62ac-4c69-8c2f-98cb0ee4503b/Tcells_TCR.rds .

# Set working directory
setwd("~/TP53_ImmuneEscape/3_DGE/")

# Delete environment variables & load favorite function
rm(list=ls())
cutf <- function(x, f=1, d="/") sapply(strsplit(x, d), function(i) paste(i[f], collapse=d))

# Load data
seu <- readRDS("~/250128_seurat_annotated_final.rds")

# Compare to other objects from late January 2025 (obsolete):
#seu_TNK <- readRDS("~/250128_Tcell_subset.rds")
#seu_TCR <- readRDS("~/Tcells_TCR.rds")
#all(colnames(seu_TNK) %in% colnames(seu))
#all(colnames(seu_TCR) %in% colnames(seu))
#df <- data.frame(row.names = levels(seu$celltype), All = rep(NA, 30), TNK = rep(NA, 30), TCR = rep(NA, 30))
#df$All <- table(seu$celltype)
#df$TNK <- table(seu_TNK$celltype)
#df$TCR <- table(seu_TCR$celltype)

# Load colors from 2.3_PvG-Colors.R
celltype_colors_df <- read.table("../celltype_colors.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE, comment.char = "")
celltype_colors <- setNames(celltype_colors_df$color, celltype_colors_df$celltype)

# Check data
seu@meta.data %>% head
seu@meta.data %>%
  sample_frac(1) %>%  # Randomly shuffle rows
      ggplot(aes(x = UMAP_1, y = UMAP_2, color = celltype)) +
      geom_scattermore(pointsize = 8, pixels = c(4096, 4096)) +
      scale_color_manual(values = celltype_colors) +
      theme_bw() +
      theme(aspect.ratio = 1,
            panel.grid = element_blank()) +
      guides(color = guide_legend(override.aes = list(size = 3)))

# Subset Seurat object to T cells
seu_T <- subset(seu, celltype %in% c("CD4 Naïve", "CD4 Effector Memory", "CD4 Memory", "Treg", "CD8 Naïve", "CD8 Effector", "CD8 Memory", "CD8 Exhausted", "γδ T"))
# Check visually
seu_T@meta.data %>% sample_frac(1) %>%
          ggplot(aes(x = UMAP_TNK_1, y = UMAP_TNK_2, color = celltype)) +
          geom_scattermore(pointsize = 16, pixels = c(4096, 4096)) +
          scale_color_manual(values = celltype_colors) +
          theme_bw() +
          theme(aspect.ratio = 1,
                panel.grid = element_blank()) +
          guides(color = guide_legend(override.aes = list(size = 3)))
ggsave("3.2.1_OurTcellAnnotationUMAP.png", width = 6, height = 3.5)

# Next, save Seurat count matrix to run TCAT (see https://github.com/immunogenomics/starCAT/blob/main/Examples/starCAT_vignette_R.ipynb)

# Extract T cell counts
counts <- LayerData(seu_T, assay = "RNA", layer = "counts")
counts[80:90, 1:10]

# Save counts matrix
dir.create("starcat")
writeMM(counts,  file = "starcat/matrix.mtx")
gzip("starcat/matrix.mtx")

# Save cell barcodes
barcodes <- colnames(counts)
write_delim(as.data.frame(barcodes),  file = "starcat/barcodes.tsv", col_names = FALSE)
gzip("starcat/barcodes.tsv")

# Save feature names
gene_names <- rownames(counts)
features <- data.frame("gene_id" = gene_names, "gene_name" = gene_names,type = "Gene Expression")
write_delim(as.data.frame(features),delim = "\t", file = "starcat/features.tsv", col_names = FALSE)
gzip("starcat/features.tsv")

# Run 3.2_PvG-TCAT.sh to generate results

# Load scores (see 3.3_PvG-TCAT_Programs.R for program analysis)
scores_tib <- read_tsv("starcat/results.scores.txt") %>% rename("cell" = "...1")

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
ggsave("3.2.2_TCAT-TcellAnnotationUMAP.png", width = 6, height = 3.5)

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
ggsave("3.2.3_Confusion_matrix.png", height = 5, width = 6)

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
ggsave("3.2.4_Side-by-side_barplot.png", height = 8, width = 5)