# Calculate the cells proportions across 2 different cohort
# Nurefsan Sariipek, 250428
# Load the libraries
library(tidyverse)
library(janitor)
library(ggrepel)
library(ggpubr)

# Empty environment
rm(list=ls())

# Set directory
setwd("~/TP53_ImmuneEscape/5_Cell_Proportions/")

# Load the seurat meta data that saved in 5.2 instead of the whole seurat object
seu_df <- read_csv("~/seu_df_250428.csv")

# Define the colors
# Celltype colors
celltype_colors_df <- read.table("~/TP53_ImmuneEscape/celltype_colors.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE, comment.char = "")
celltype_colors <- setNames(celltype_colors_df$color, celltype_colors_df$celltype)

# Cohort colors
cohort_colors <- c("long-term-remission" = "#546fb5FF","relapse" = "#e54c35ff")

# Levels disapear after turning seu object to metadata, add the new levels
my_levels <- c("HSC MPP", "MEP", "EoBasoMast Precursor", "Megakaryocyte Precursor", "Early Erythroid", "Late Erythroid", 
               "LMPP", "Cycling Progenitor", "Early GMP", "Late GMP", "Pro-Monocyte", "CD14 Mono", 
               "CD16 Mono", "cDC", "pDC", "Early Lymphoid", "Pro-B", "Pre-B", 
               "Immature B", "Mature B", "Plasma Cell", "CD4 Naive", "CD4 Central Memory", "CD4 Effector Memory", 
               "CD4 Regulatory", "CD8 Naive", "CD8 Central Memory", "CD8 Effector Memory 1", "CD8 Effector Memory 2", 
               "CD8 Tissue Resident Memory", "T Proliferating", "NK", "NK CD56high", "NK Proliferating", "Stromal")

# Turn into right format for plotting
bar_data <- seu_df %>%
  filter(library_type == "MNC",timepoint %in% c("3","5","6") , sample_status =="remission") %>% 
  drop_na(celltype) %>%
  group_by(patient_id, celltype) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(patient_id) %>%
  mutate(percent = count / sum(count) * 100)

bar_data$celltype <- factor(bar_data$celltype, levels = my_levels)

# Plot
p1 <- ggplot(bar_data, aes(x = patient_id, y = percent, fill = celltype)) +
  geom_bar(stat = "identity", width = 0.8) +
  labs(x = "Patient ID", y = "Cell Type Proportion (%)") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = celltype_colors) +
  theme_minimal(base_size = 9) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8, color = "black"),
    axis.text.y = element_text(size = 8, color = "black"),
    axis.title = element_text(size = 8, color = "black"),
    legend.title = element_blank(),
    panel.grid= element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"))
p1

# Save as a pdf
pdf("5.1_Allcells_proportions_3-6mo.pdf", width = 10, height = 5)
p1
dev.off()

##################################################################################################

# Calculate total cells and myeloid cells per patient ==percentage of all myeloid cells out of all cell types in each sample (MNC libraries only)
proportions_df <- seu_df %>%
  filter(timepoint %in% c("3", "5", "6"),
         sample_status == "remission", library_type=="MNC") %>%
  drop_na(celltype) %>%
  group_by(patient_id, cohort) %>%
  mutate(total_cells = n()) %>%  # total cells per sample
  count(patient_id, cohort, total_cells, celltype, name = "cell_count") %>%
  mutate(percent = (cell_count / total_cells) * 100) 

proportions_df$celltype <- factor(proportions_df$celltype, levels = my_levels)

# Plot 
p2 <- ggplot(proportions_df, aes(x = cohort, y = percent, fill = cohort)) +
  geom_boxplot(width = 0.7, outlier.shape = NA, alpha = 0.9) +
  geom_jitter(shape = 21, size = 2, stroke = 0.2, color = "black", width = 0.15) +
  facet_wrap(~ celltype, ncol = 9) +
  scale_fill_manual(values = cohort_colors) +
  labs(y = "% of all cells", x = NULL) +
  ggpubr::stat_compare_means(aes(x = cohort, y = percent, group = cohort),
    method = "t.test",
    label = "p.format",
    hide.ns = TRUE,
    label.y = 30) +
  theme_minimal(base_size = 6) +
  theme(axis.text.x = element_text(size = 8, color = "black", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 8, color = "black"),
    axis.title.y = element_text(size =11, face = "plain", color = "black"),
    strip.text = element_text(size = 8, face = "plain", color = "black"),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
    legend.position = "none",
    aspect.ratio = 2,
    plot.margin = margin(4, 4, 4, 4))

p2
# Save as a pdf
pdf("5.1_Allcells_per_celltype_proportions_3-6mo_pvalue.pdf", width = 10, height = 8)
p2
dev.off()

########################

