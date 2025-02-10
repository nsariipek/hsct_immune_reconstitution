# Peter van Galen, 250206
# Assess TCAT programs in T cells

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

seu@meta.data %>% filter(UMAP_TNK_1 != "NA", UMAP_TNK_2 != "NA") %>%
      sample_frac(1) %>%  # Randomly shuffle rows
          ggplot(aes(x = UMAP_TNK_1, y = UMAP_TNK_2, color = celltype)) +
          geom_scattermore(pointsize = 8, pixels = c(4096, 4096)) +
          scale_color_manual(values = celltype_colors) +
          theme_bw() +
          theme(aspect.ratio = 1,
                panel.grid = element_blank()) +
          guides(color = guide_legend(override.aes = list(size = 3)))

# Now, using part of the StarCAT tutorial, save output files to run TCAT (see https://github.com/immunogenomics/starCAT/blob/main/Examples/starCAT_vignette_R.ipynb)

# Subset Seurat object to T cells
seu_T <- subset(seu, celltype %in% c("CD4 Naïve", "CD4 Effector Memory", "CD4 Memory", "Treg", "CD8 Naïve", "CD8 Effector", "CD8 Memory", "CD8 Exhausted", "γδ T"))
# Check visually
seu_T@meta.data %>% sample_frac(1) %>%
          ggplot(aes(x = UMAP_TNK_1, y = UMAP_TNK_2, color = celltype)) +
          geom_scattermore(pointsize = 8, pixels = c(4096, 4096)) +
          scale_color_manual(values = celltype_colors) +
          theme_bw() +
          theme(aspect.ratio = 1,
                panel.grid = element_blank()) +
          guides(color = guide_legend(override.aes = list(size = 3)))
# Extract T cell counts
counts <- LayerData(seu_T, assay = "RNA", layer = "counts")
counts[80:90, 1:10]

# Output counts matrix
dir.create("starcat")
writeMM(counts,  file = "starcat/matrix.mtx")
gzip("starcat/matrix.mtx")

# Output cell barcodes
barcodes <- colnames(counts)
write_delim(as.data.frame(barcodes),  file = "starcat/barcodes.tsv", col_names = FALSE)
gzip("starcat/barcodes.tsv")

# Output feature names
gene_names <- rownames(counts)
features <- data.frame("gene_id" = gene_names, "gene_name" = gene_names,type = "Gene Expression")
write_delim(as.data.frame(features),delim = "\t", file = "starcat/features.tsv", col_names = FALSE)
gzip("starcat/features.tsv")

# Run 3.2_PvG-TCAT.sh to generate results

# Load results
usage_tib <- read_tsv("starcat/results.rf_usage_normalized.txt") %>% rename("cell" = "...1")
scores_tib <- read_tsv("starcat/results.scores.txt") %>% rename("cell" = "...1")

# Compare cell type annotations with Multinomial_Label from scores
metadata_tib <- as_tibble(seu_T@meta.data, rownames = "cell")
metadata_tib <- left_join(metadata_tib, scores_tib)

# UMAP
metadata_tib %>% sample_frac(1) %>%
          ggplot(aes(x = UMAP_TNK_1, y = UMAP_TNK_2, color = Multinomial_Label)) +
          geom_scattermore(pointsize = 8, pixels = c(4096, 4096)) +
          theme_bw() +
          theme(aspect.ratio = 1,
                panel.grid = element_blank()) +
          guides(color = guide_legend(override.aes = list(size = 3)))

# Heatmap
metadata_tib %>% select(celltype, Multinomial_Label) %>%
      group_by(celltype, Multinomial_Label) %>% count() %>%
      ggplot(aes(x = Multinomial_Label, y = celltype, fill = n)) +
      geom_tile() +
      scale_fill_gradient(low = "white", high = "red") + 
      theme_minimal() +
      labs(x = "Multinomial Label", y = "Cell Type", fill = "Count") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Evaulate program usages
metadata_tib <- left_join(metadata_tib, usage_tib)

# Make a list of relevant programs
programs <- setdiff(colnames(usage_tib), "cell")
programs <- setdiff(programs, programs[grepl("Doublet", programs)])
programs <- setdiff(programs, "Poor-Quality", )
# Remove programs that resemble cell types (I'm not sure there is a good rationale for this)
programs <- setdiff(programs, unique(gsub("_", "-", scores$Multinomial_Label)))

# Filter data for relevant time points etc.
metadata_subset_tib <- metadata_tib %>%
      filter(timepoint %in% c("3","5","6"),  # 100 days after transplant
             sample_status == "remission") %>%
      mutate(cohort_binary = substr(cohort, 3, nchar(cohort))) %>%
      select(patient_id, celltype, cohort_binary, programs)

# Calculate average program scores per patient
metadata_subset_tib %>%
      pivot_longer(cols = programs, names_to = "program") %>% # Footnote 1
      group_by(patient_id, cohort_binary, celltype, program) %>%
      summarize(mean_pp = mean(value), .groups = "drop") %>% # Mean program usage per patient (2)
      group_by(cohort_binary, celltype, program) %>%
      summarize(mean_p = mean(mean_pp), .groups = "drop") %>% # Mean program usage per cohort (3)
      ggplot(aes(x = celltype, y = program, fill = mean_p)) +
      geom_tile() +
      scale_fill_gradient(low = "white", high = "red") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      facet_wrap(~ cohort_binary)
# 1. At this point, the tibble has 118,503 cells x 52 programs = 6,162,156 rows
# 2. At this point, the tibble has 28 patients x 9 cell types x 52 programs = 13,104 rows
# 3. At this point, the tible has 2 cohorts x 9 cell types x 52 programs = 936 rows

ggsave("3.2_Program-usage-cohorts.png", width = 8, height = 8)

# That looks pretty similar. Still: look at differences
non_relapsed_tib <- metadata_subset_tib %>% filter(cohort_binary == "Non-relapsed") %>%
      pivot_longer(cols = programs, names_to = "program") %>% # Footnote 1
      group_by(patient_id, cohort_binary, celltype, program) %>%
      summarize(mean_pp = mean(value), .groups = "drop") %>%
      group_by(celltype, program) %>%
      summarize(mean_p = mean(mean_pp), .groups = "drop")

relapsed_tib <- metadata_subset_tib %>% filter(cohort_binary == "Relapsed") %>%
      pivot_longer(cols = programs, names_to = "program") %>% # Footnote 1
      group_by(patient_id, cohort_binary, celltype, program) %>%
      summarize(mean_pp = mean(value), .groups = "drop") %>%
      group_by(celltype, program) %>%
      summarize(mean_p = mean(mean_pp), .groups = "drop")

join_tib <- full_join(non_relapsed_tib, relapsed_tib, by = c("celltype", "program"),
                      suffix = c(".non-relapsed", ".relapsed"))
join_tib <- join_tib %>% mutate(diff = mean_p.relapsed - `mean_p.non-relapsed`)

join_tib %>%
      ggplot(aes(x = celltype, y = program, fill = diff)) +
      geom_tile() +
      scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 5, name = "RdBu"))) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("3.2_Program-usage-diff.png", width = 8, height = 12)

