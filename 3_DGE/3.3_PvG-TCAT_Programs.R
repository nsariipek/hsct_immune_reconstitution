# Peter van Galen, 250206
# Assess TCAT programs in T cells

library(Seurat)
library(tidyverse)
library(scattermore)

# Set working directory
setwd("~/TP53_ImmuneEscape/3_DGE/")

# Delete environment variables & load favorite function
rm(list=ls())
cutf <- function(x, f=1, d="/") sapply(strsplit(x, d), function(i) paste(i[f], collapse=d))

# Load data
seu <- readRDS("~/250128_seurat_annotated_final.rds")

# Load colors from 2.3_PvG-Colors.R
celltype_colors_df <- read.table("../celltype_colors.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE, comment.char = "")
celltype_colors <- setNames(celltype_colors_df$color, celltype_colors_df$celltype)

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

# Combine TCAT program usage results with Seurat metadata
usage_tib <- read_tsv("starcat/results.rf_usage_normalized.txt") %>% rename("cell" = "...1")
metadata_tib <- as_tibble(seu_T@meta.data, rownames = "cell")
metadata_tib <- left_join(metadata_tib, usage_tib)

# Make a list of activity programs from the paper (Figure 2C)
programs <- c("BCL2/FAM13A", "TIMD4/TIM3", "Th2-Activated", "ICOS/CD38", "Th17-Activated", "CD172a/MERTK", "SOX4/TOX2", "Heatshock", "OX40/EBI3", 
              "CD40LG/TXNIP", "Exhaustion", "CellCycle-G2M", "CellCycle-S", "Cytoskeleton", "Cytotoxic", "ISG", "CellCycle-Late-S", "HLA", 
              "NME1/FABP5", "Translation", "IEG", "IL10/IL19", "RGCC/MYADM", "Metallothionein", "Multi-Cytokine", "CTLA4/CD38")

# Filter data for relevant time points etc.
metadata_subset_tib <- metadata_tib %>%
      filter(timepoint %in% c("3","5","6"),  # 100 days after transplant
             sample_status == "remission") %>%
      mutate(cohort_binary = substr(cohort, 3, nchar(cohort))) %>%
      select(patient_id, celltype, cohort_binary, all_of(programs))

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

ggsave("3.3.1_Program-usage-cohorts.png", width = 8, height = 8)

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

ggsave("3.3.2_Program-usage-diff.png", width = 5, height = 8)

