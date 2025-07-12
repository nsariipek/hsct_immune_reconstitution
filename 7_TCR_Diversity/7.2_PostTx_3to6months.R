# Nurefsan Sariipek and Peter van GAlen
# Date: January 22nd, 2024
# Updated: June 6, 2025
# Analyze post 3-6 months remission samples and subset donor/host CD4/CD8 compartments using subsampling based on cell numbers which is different from scRepertoire built-in subsampling

# Load the libraries
library(tidyverse)
library(Seurat)
library(ggpubr)

# Empty environment
rm(list=ls())

# Set working directory for Nurefsan/Terra:
setwd("~/TP53_ImmuneEscape/7_TCR_Diversity/")
# For Peter:
setwd("~/DropboxMGB/Projects/ImmuneEscapeTP53/TP53_ImmuneEscape/7_TCR_Diversity")

# Load necessary functions
cutf <- function(x, f=1, d="/") sapply(strsplit(x, d), function(i) paste(i[f], collapse=d))
source("DiversityFunctions.R")

# Load final Seurat object including TCR calls
seu <- readRDS("../AuxiliaryFiles/250528_Seurat_complete.rds")


######## PREPARE DATA ########

# Keep only annotated T cell clusters
TCAT_cells <- colnames(seu)[!is.na(seu$TCAT_Multinomial_Label)]
seu_T <- subset(seu, cells = TCAT_cells)

# Turn metadata to a tibble and keep only needed variables
metadata_tib <- as_tibble(seu_T@meta.data, rownames = "cell")
metasubset_tib <- metadata_tib %>% select(cell, orig.ident, cohort, patient_id, timepoint, sample_status, TP53_status, TCAT_Multinomial_Label, CTstrict, souporcell_origin)

# Keep only 3-6M remission samples and exclude early relapse cohort
metasubset2_tib <- metasubset_tib %>% filter(between(timepoint, 3, 6),
  sample_status == "remission",
  ! patient_id %in% c("P30", "P31", "P32", "P33"))

# Add group for TCR diversity calculation
metasubset2_tib$group <- paste0(metasubset2_tib$patient_id, "_", metasubset2_tib$sample_status, "_3-6M")
# Check
metasubset2_tib %>% group_by(patient_id, sample_status, timepoint, group) %>% count() # %>% view


######## SUBSET ########

# Subsetting for different cell types. Choose one!

# All cells
cell_subset <- "_all"
metasubset3_tib <- metasubset2_tib

# CD4+ T cells
cell_subset <- "_CD4"
metasubset3_tib <- metasubset2_tib %>% subset(TCAT_Multinomial_Label %in% c("CD4_Naive", "CD4_CM", "CD4_EM", "Treg"))

# CD8+ T cells
cell_subset <- "_CD8"
metasubset3_tib <- metasubset2_tib %>% subset(TCAT_Multinomial_Label %in% c("CD8_Naive", "CD8_CM", "CD8_EM", "CD8_TEMRA"))

# Recipient cells
cell_subset <- "_recipient"
metasubset3_tib <- metasubset2_tib %>% subset(souporcell_origin == "recipient")

# Donor cells
cell_subset <- "_donor"
metasubset3_tib <- metasubset2_tib %>% subset(souporcell_origin == "donor")


######## CALCULATE AND PLOT DIVERSITY ########

# Cell threshold
min_cells <- 300

# Subset to samples with a reasonable number of cells
metasubset3_tib %>% group_by(group) %>% count()
metasubset3_tib %>% group_by(group) %>% count() %>%
  ggplot(aes(x = group, y = n)) +
  geom_point() +
  geom_hline(yintercept = min_cells, color = "red") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank())

# Filter for samples with sufficient cells
metasubset4_tib <- metasubset3_tib %>% group_by(group) %>% filter(n() >= min_cells)

# Sample min_cells from each group. This gives all the plots below the same y-axis range
metasubset4_tib <- metasubset4_tib %>% group_by(group) %>% slice_sample(n = min_cells)

# Split data into list
metasubset_ls <- split(metasubset4_tib, f = metasubset4_tib$group)

# Calculate diversity (works only with lists)
diversities_df = compute_diversity(metasubset_ls, "CTstrict", 1000)
diversities_df
max(diversities_df$inv.simpson)

# Add information to make annotated plots
metasubset_summary <- metasubset4_tib %>% group_by(group, cohort, patient_id, TP53_status) %>% summarize
joined_tibble <- diversities_df %>% left_join(metasubset_summary)

# Visualize barplot
plot_TCR <- function(df, title = NULL) {
  df %>%
  ggplot(aes(x = cohort, y = inv.simpson, fill = cohort)) +
  geom_bar(stat = "summary", fun = mean, width = 0.6, color = "black") +
  geom_jitter(width = 0.25, size = 4, alpha = 0.7, color = "black") +
  stat_summary(fun.data = mean_se, geom = "errorbar", 
               width = 0.2, linewidth = 1, color = "black") +
  scale_fill_manual(values = c("long-term-remission" = "#546fb5FF","relapse" = "#e54c35ff")) +
  coord_cartesian(ylim = c(0,300)) +
  stat_compare_means(
    aes(group = cohort), 
    method = "wilcox.test", 
    label = "p.format",
    label.x.npc = "center",
    label.x = 1.5,
    label.y = 285,
    size = 5) +
    theme_pubr(base_size = 16) +
  labs(y = "Inverse Simpson Index",
       x = "Cohort") +
  theme(aspect.ratio = 2,
    axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 14),
    axis.title.y = element_text(size = 16),
    axis.title.x = element_text(size = 16),
    legend.position = "none",
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12))}

TCR_all <- plot_TCR(joined_tibble, title = "All samples")
TCR_MT <- plot_TCR(filter(joined_tibble, TP53_status == "MUT"), title = "TP53-MUT samples")
TCR_WT <- plot_TCR(filter(joined_tibble, TP53_status == "WT"), title = "TP53-WT samples") 

# Save as a pdf: combined TP53-MUT and TP53-WT, or separate
pdf(paste0("7.2_Post-transplant", cell_subset, "_both.pdf"), width = 4, height = 6)
TCR_all
dev.off()

pdf(paste0("7.2_Post-transplant", cell_subset, "_MT.pdf"), width = 4, height = 6)
TCR_MT
dev.off()

pdf(paste0("7.2_Post-transplant", cell_subset, "_WT.pdf"), width = 4, height = 6)
TCR_WT
dev.off()

