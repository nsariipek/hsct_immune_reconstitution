# Peter van Galen, 250408, updated 250504
# Run TCAT and evaluate the algorithm's cell type labels

# Load libraries
library(Seurat)
library(tidyverse)
library(scattermore)
library(Matrix)
library(R.utils)

# Set working directory
repo_root <- system("git rev-parse --show-toplevel", intern = TRUE)
setwd(paste0(repo_root, "/05_DGE"))

# Clear environment variables & load favorite function
rm(list = ls())
cutf <- function(x, f = 1, d = "/") {
      sapply(strsplit(x, d), function(i) paste(i[f], collapse = d))
}

# Load data
seu <- readRDS("../AuxiliaryFiles/250426_Seurat_annotated.rds")

# Subset for T cells
T_celltypes <- c(
      "CD4 Naive",
      "CD4 Central Memory",
      "CD4 Effector Memory",
      "CD4 Regulatory",
      "CD8 Naive",
      "CD8 Central Memory",
      "CD8 Effector Memory 1",
      "CD8 Effector Memory 2",
      "CD8 Tissue Resident Memory",
      "T Proliferating"
)
seu_T <- subset(seu, subset = celltype %in% T_celltypes)

# Load colors
celltype_colors_df <- read.table(
      "../celltype_colors.txt",
      sep = "\t",
      header = T,
      stringsAsFactors = F,
      comment.char = ""
)
celltype_colors <- setNames(
      celltype_colors_df$color,
      celltype_colors_df$celltype
)

# Check data visually
DimPlot(seu_T, group.by = "celltype") +
      scale_color_manual(values = celltype_colors) +
      coord_cartesian(xlim = c(-13, -5), ylim = c(-7, 7))

# Next, save Seurat count matrix to run TCAT (see https://github.com/immunogenomics/starCAT/blob/main/Examples/starCAT_vignette_R.ipynb)

# Extract T cell counts
counts <- LayerData(seu_T, assay = "RNA", layer = "counts")
counts[80:90, 1:10]

# Save counts matrix
dir.create("../AuxiliaryFiles/starCAT")
writeMM(counts, file = "../AuxiliaryFiles/starCAT/matrix.mtx")
gzip("../AuxiliaryFiles/starCAT/matrix.mtx", overwrite = T, remove = T)

# Save cell barcodes
barcodes <- colnames(counts)
write_delim(
      as.data.frame(barcodes),
      file = "../AuxiliaryFiles/starCAT/barcodes.tsv",
      col_names = FALSE
)
gzip("../AuxiliaryFiles/starCAT/barcodes.tsv", overwrite = T, remove = T)

# Save feature names
gene_names <- rownames(counts)
features <- data.frame(
      "gene_id" = gene_names,
      "gene_name" = gene_names,
      type = "Gene Expression"
)
write_delim(
      as.data.frame(features),
      delim = "\t",
      file = "../AuxiliaryFiles/starCAT/features.tsv",
      col_names = FALSE
)
gzip("../AuxiliaryFiles/starCAT/features.tsv", overwrite = T, remove = T)

# Run 3.1_PvG-TCAT.sh to generate results (this was done on a Google Cloud Platform virtual machine terminal). This will create a folder 5.1_starCAT with the results
#cd /home/unix/vangalen/hsct_immune_reconstitution/05_DGE
#mkdir 5.1_starCAT
#./5.1_starCAT.sh
#rm -r cache

# Zip the results to enable syncing over GitHub
gzip("5.1_starCAT/starCAT.rf_usage_normalized.txt", overwrite = T, remove = T)
gzip("5.1_starCAT/starCAT.scores.txt", overwrite = T, remove = T)

# The rest of the script compares TCAT cell annotations to our own and is not used for the paper. Start with loading scores.
scores_tib <- read_tsv("5.1_starCAT/starCAT.scores.txt.gz") %>%
      rename("cell" = "...1")

# Compare cell type annotations with Multinomial_Label from scores
metadata_tib <- as_tibble(seu_T@meta.data, rownames = "cell")
metadata_tib$UMAP_1 <- seu_T@reductions$umap_bmm@cell.embeddings[, 1]
metadata_tib$UMAP_2 <- seu_T@reductions$umap_bmm@cell.embeddings[, 2]
metadata_tib <- left_join(metadata_tib, scores_tib)

# Wrangle and assign colors
metadata_tib$Multinomial_Label <- factor(
      metadata_tib$Multinomial_Label,
      levels = c(
            "CD4_Naive",
            "CD4_CM",
            "CD4_EM",
            "Treg",
            "CD8_Naive",
            "CD8_CM",
            "CD8_EM",
            "CD8_TEMRA",
            "gdT",
            "MAIT"
      )
)

# UMAP
metadata_tib %>%
      sample_frac(1) %>%
      ggplot(aes(x = UMAP_1, y = UMAP_2, color = Multinomial_Label)) +
      geom_scattermore(pointsize = 16, pixels = c(4096, 4096)) +
      coord_cartesian(xlim = c(-13, -5), ylim = c(-7, 7)) +
      scale_color_manual(values = celltype_colors) +
      theme_bw() +
      theme(aspect.ratio = 1, panel.grid = element_blank()) +
      guides(color = guide_legend(override.aes = list(size = 3)))

# Compare cell type numbers
p1 <- metadata_tib %>%
      count(celltype) %>%
      ggplot(aes(x = celltype, y = n, fill = celltype)) +
      geom_bar(stat = "identity") +
      scale_fill_manual(values = celltype_colors) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
p2 <- metadata_tib %>%
      count(Multinomial_Label) %>%
      ggplot(aes(x = Multinomial_Label, y = n, fill = Multinomial_Label)) +
      geom_bar(stat = "identity") +
      scale_fill_manual(values = celltype_colors) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
p1 + p2

# Heatmap
df <- data.frame(
      celltype = rep(1:10, each = 10),
      Multinomial_Label = rep(1:10, times = 10)
)
metadata_tib %>%
      select(celltype, Multinomial_Label) %>%
      count(celltype, Multinomial_Label) %>%
      complete(celltype, Multinomial_Label, fill = list(n = 0)) %>%
      #group_by(celltype) %>%
      #mutate(prop = n / sum(n)) %>%
      ggplot(aes(x = Multinomial_Label, y = celltype, fill = n)) +
      geom_tile(color = "grey", size = 0.25) +
      geom_tile(
            data = df[df$celltype == df$Multinomial_Label, ],
            color = "black",
            size = 0.5,
            fill = NA
      ) +
      scale_fill_gradient(low = "white", high = "red") +
      theme_minimal() +
      labs(x = "Multinomial Label", y = "Cell Type", fill = "Count") +
      theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
            aspect.ratio = 10 / 10,
            panel.grid = element_blank()
      )
ggsave("3.1_TCAT_confusion.pdf")

# Compare annotations side by side (stacked cells)
metadata_tib %>%
      select(cell, celltype, Multinomial_Label) %>%
      arrange(desc(celltype), desc(Multinomial_Label)) %>%
      pivot_longer(
            cols = -cell,
            names_to = "annotation_approach",
            values_to = "anno"
      ) %>%
      mutate(cell = factor(cell, levels = unique(cell))) %>%
      ggplot(aes(x = annotation_approach, y = cell, fill = anno)) +
      geom_tile() +
      scale_fill_manual(values = celltype_colors) +
      theme_minimal() +
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

# Compare CD4/CD8 ratio (similar to 3.3_CD4-CD8_ratio.R but with TCAT labels)
metadata_tib %>%
      filter(
            sample_status == "remission",
            timepoint %in% c(3, 5, 6),
            TP53_status == "MUT"
      ) %>%
      mutate(
            type = case_when(
                  grepl("CD4|Treg", Multinomial_Label) ~ "CD4",
                  grepl("CD8", Multinomial_Label) ~ "CD8"
            )
      ) %>%
      # group_by(Multinomial_Label, type) %>% count # check
      group_by(patient_id, cohort, type) %>%
      count() %>%
      ungroup() %>%
      group_by(patient_id) %>%
      pivot_wider(names_from = type, values_from = n) %>%
      mutate(ratio = CD4 / CD8) %>%
      ggplot(aes(x = cohort, y = ratio, fill = cohort)) +
      geom_boxplot(width = 0.7, alpha = 0.9, outlier.shape = NA) +
      geom_jitter(shape = 21, size = 2, color = "black") +
      scale_fill_manual(
            values = c(
                  "long-term-remission" = "#546fb5",
                  "relapse" = "#e54c35"
            )
      ) +
      labs(y = "CD4/CD8 T cell ratio", x = NULL) +
      stat_compare_means(
            aes(group = cohort),
            method = "wilcox.test",
            label = "p.format",
            label.y = 4,
            label.x = 1.2,
            size = 3
      ) +
      theme_minimal() +
      theme(
            aspect.ratio = 2,
            panel.grid = element_blank(),
            axis.text = element_text(color = "black"),
            axis.text.x = element_text(
                  color = "black",
                  angle = 45,
                  hjust = 1
            ),
            panel.border = element_rect(color = "black", fill = NA),
            axis.ticks = element_line(color = "black"),
            legend.position = "none"
      )
