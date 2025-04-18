# Nurefsan Sariipek, 231103 updated at 250416
# Check T/NK cell proportions
# Load the libraries
library(tidyverse)
library(janitor)
library(RColorBrewer)
library(randomcoloR)
library(ggrepel)
library(ggpubr)
library(ggnewscale)

# Empty environment
rm(list=ls())

# For Nurefsan:
setwd("~/TP53_ImmuneEscape/5_Cell_Proportions/")

# # Read the latest and turn into data frame
# seu <- readRDS("~/250418_Seurat_all_cells_annotated.rds")
# seu_df <- seu@meta.data
# #Save as this 
# write_csv(seu_df,"~/seu_df_250418.csv")

# Load the saved seu dataframe
seu_df <- read_csv("~/seu_df_250418.csv")

# Define the color palette
# Celltype colors
celltype_colors_df <- read.table("~/TP53_ImmuneEscape/celltype_colors.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE, comment.char = "")
celltype_colors <- setNames(celltype_colors_df$color, celltype_colors_df$celltype)

# Cohort colors
cohort_colors <- c("long-term-remission" = "#546fb5FF","relapse" = "#e54c35ff")

# Levels disapear after turning seu object to metadata, add the new levels
my_levels <- c("CD4 Naive", "CD4 Memory", "CD4 Effector Memory", "Treg", 
               "CD8 Naive", "CD8 Memory", "CD8 Effector", "CD8 Exhausted", "Gamma-Delta T", 
               "NK-T", "Adaptive NK", "CD56 Bright NK", "CD56 Dim NK", "Cycling T-NK")


# Prepare the data
bar_data <- seu_df %>%
mutate(patient_id = factor(patient_id)) %>%
filter(timepoint %in% c("3","5","6") , 
       sample_status =="remission", 
       celltype %in% c("CD4 Naive","CD4 Memory","CD4 Effector Memory","Treg","CD8 Naive","CD8 Memory", "CD8 Effector","CD8 Exhausted", "Gamma-Delta T")) %>%
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
    panel.grid.major.x = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"))
p1

# Save as a pdf
pdf("5.2_Tcells_proportions_3-6mo.pdf", width = 10, height = 5)
p1
dev.off()

###############################################################
# Calculate the proportions for T cells in MUT samples

proportions_df <- seu_df %>%
  filter(celltype %in% c("CD4 Naive", "CD4 Memory", "CD4 Effector Memory", "Treg",
                         "CD8 Naive", "CD8 Memory", "CD8 Effector", "CD8 Exhausted",
                         "Gamma-Delta T"),
         timepoint %in% c("3", "5", "6"),
         sample_status == "remission",
         TP53_status == "MUT") %>%
  group_by(patient_id, cohort) %>%
  count(celltype, name = "n") %>%   # Simpler and more reliable than `tabyl()`
  mutate(total_T_cells = sum(n),
         percent_within_T = (n / total_T_cells) * 100) %>%
  ungroup()

proportions_df$celltype <- factor(proportions_df$celltype, levels = my_levels)

# Plot
p2 <- ggplot(proportions_df, aes(x = cohort, y = percent_within_T, fill = cohort)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.9, linewidth = 0.5, color = "black") + 
  geom_jitter(shape = 21, size = 1.8, width = 0.15, stroke = 0.2, color = "black", aes(fill = cohort)) +
  coord_cartesian(ylim = c(0, 50)) +
  scale_y_continuous(breaks = seq(0, 50, 10)) +  # Add y-axis ticks every 10
  facet_wrap(~ celltype, ncol = 10) +
  scale_fill_manual(values = cohort_colors) +
  labs(y = "% within total T cells", x = NULL) +
  stat_compare_means(
    aes(x = cohort, y = percent_within_T, group = cohort),
    method = "wilcox.test",
    label = "p.format",
    size = 4,
    hide.ns = TRUE) +
  theme_minimal(base_size = 8) +
  theme(
    axis.text.x = element_text(size = 8, color = "black", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 8, color = "black"),
    axis.title.y = element_text(size = 8, face = "plain", color = "black"),
    strip.text = element_text(size = 8, face = "plain", color = "black"),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
    axis.ticks = element_line(color = "black", linewidth = 0.3),
    axis.ticks.length = unit(2, "pt"),
    legend.position = "none",
    aspect.ratio = 2,
    plot.margin = margin(4, 4, 4, 4))

# Check the plot
p2

# Save as a pdf
pdf("5.2_T_proportions_TP53_MT_pvalue.pdf", width = 8, height = 6)
p2
dev.off()

# Calculate the proportions for T cells in WT samples
proportions_df <- seu_df %>%
  filter(celltype %in% c("CD4 Naive", "CD4 Memory", "CD4 Effector Memory", "Treg",
                         "CD8 Naive", "CD8 Memory", "CD8 Effector", "CD8 Exhausted",
                         "Gamma-Delta T"),
         timepoint %in% c("3", "5", "6"),
         sample_status == "remission",
         TP53_status == "WT") %>%
  group_by(patient_id, cohort) %>%
  count(celltype, name = "n") %>%   # Simpler and more reliable than `tabyl()`
  mutate(total_T_cells = sum(n),
         percent_within_T = (n / total_T_cells) * 100) %>%
  ungroup()

proportions_df$celltype <- factor(proportions_df$celltype, levels = my_levels)


# Plot
p3 <- ggplot(proportions_df, aes(x = cohort, y = percent_within_T, fill = cohort)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.9, linewidth = 0.5, color = "black") + 
  geom_jitter(shape = 21, size = 1.8, width = 0.15, stroke = 0.2, color = "black", aes(fill = cohort)) +
  coord_cartesian(ylim = c(0, 50)) +
  scale_y_continuous(breaks = seq(0, 50, 10)) +  # Add y-axis ticks every 10
  facet_wrap(~ celltype, ncol = 10) +
  scale_fill_manual(values = cohort_colors) +
  labs(y = "% within total T cells", x = NULL) +
  stat_compare_means(
    aes(x = cohort, y = percent_within_T, group = cohort),
    method = "wilcox.test",
    label = "p.format",
    size = 4,
    hide.ns = TRUE) +
  theme_minimal(base_size = 8) +
  theme(
    axis.text.x = element_text(size = 8, color = "black", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 8, color = "black"),
    axis.title.y = element_text(size = 8, face = "plain", color = "black"),
    strip.text = element_text(size = 8, face = "plain", color = "black"),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
    axis.ticks = element_line(color = "black", linewidth = 0.3),
    axis.ticks.length = unit(2, "pt"),
    legend.position = "none",
    aspect.ratio = 2,
    plot.margin = margin(4, 4, 4, 4))

# Check the plot
p3

# Save as a pdf
pdf("5.2_T_proportions_TP53_WT_pvalue.pdf", width = 8, height = 6)
p3
dev.off()




