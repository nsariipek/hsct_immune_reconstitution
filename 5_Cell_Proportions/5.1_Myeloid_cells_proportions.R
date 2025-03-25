# Calculate the myeloid cells proportions
# Load the libraries
library(tidyverse)
library(janitor)
library(ggrepel)
library(ggpubr)

# Empty environment
rm(list=ls())

# set directory
setwd("~/TP53_ImmuneEscape/5_Cell_Proportions/")

# Load the saved dataframe that contains the information about souporcell information
final_df <- read_csv("~/final_dataset.csv")

# Calculate the proportions of each cell type within the myeloid lineage(including MIX+MNC libraries)
  proportions_df1 <- final_df %>%
  filter(celltype %in% c("Progenitors","Early Erythroids","Mid Erythroids" ,"Late Erythroids","cDC","Pro Monocytes", "Monocytes","Non Classical Monocytes") &
           timepoint %in% c("3","5","6") & sample_status =="remission" 
          # origin == "donor"
  ) %>%
  group_by(sample_id, survival) %>% reframe(tabyl(celltype)) %>%
  mutate(percent = percent*100) %>%
  mutate(celltype = factor(celltype,
                           levels = c("Progenitors","Early Erythroids","Mid Erythroids" ,"Late Erythroids","pDC","cDC","Pro Monocytes", "Monocytes","Non Classical Monocytes")))

# Visualize 
p1 <- ggplot(proportions_df1, aes(x = survival, y = percent, fill = survival)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.9, linewidth = 0.5, color = "black") +  # Thicker lines
  geom_jitter(shape = 21, size = 1.8, width = 0.15, stroke = 0.2, color = "black", aes(fill = survival)) +
  coord_cartesian(ylim = c(0, 60)) +
  facet_wrap(~ celltype, ncol = 4) +  
  scale_fill_manual(values = c("Relapsed" = "#E64B35FF", "Non-relapsed" = "#4DBBD5FF")) +
  labs(y = "% of cells", x = NULL) +
  stat_compare_means(
    method = "wilcox.test",
    label = "p.format",
    size = 2.3,
    hide.ns = FALSE) +
  theme_minimal(base_size = 8) +
  theme(
    axis.text.x = element_text(size = 12, color = "black", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 11, color = "black"),
    axis.title.y = element_text(size = 15, face = "plain"),
    strip.text = element_text(size = 10, face = "bold"),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
    legend.position = "none",
    aspect.ratio = 2,  # ← Keeps it narrow/tall
    plot.margin = margin(4, 4, 4, 4))
p1

# Save as a pdf
pdf("5.1_Myeloid_proportions_3-6mo.pdf", width = 10, height = 8)
p1
dev.off()

######################################################################################################################


# Calculate total cells and myeloid cells per sample ==  percentage of all myeloid cells out of all cell types in each sample 

myeloid_types <- c("Progenitors","Early Erythroids","Mid Erythroids",
                   "Late Erythroids","pDC","Pro Monocytes", 
                   "Monocytes","Non Classical Monocytes")

proportions_df2 <- final_df %>%
  filter(timepoint %in% c("3", "5", "6"),
         sample_status == "remission", library_type=="MNC") %>%
  group_by(sample_id, survival) %>%
  mutate(total_cells = n()) %>%  # total cells per sample
  filter(celltype %in% myeloid_types) %>%
  count(sample_id, survival, total_cells, celltype, name = "cell_count") %>%
  mutate(percent = (cell_count / total_cells) * 100) %>%
  mutate(celltype = factor(celltype,
                           levels = myeloid_types))

p2 <- ggplot(proportions_df2, aes(x = survival, y = percent, fill = survival)) +
  geom_boxplot(width = 0.7, outlier.shape = NA, alpha = 0.9) +
  geom_jitter(shape = 21, size = 2, stroke = 0.2, color = "black", width = 0.15) +
  facet_wrap(~ celltype, ncol = 4) +
  scale_fill_manual(values = c("Relapsed" = "#E64B35FF", "Non-relapsed" = "#4DBBD5FF")) +
  labs(y = "% of all cells", x = NULL) +
  stat_compare_means(method = "wilcox.test", label = "p.signif", hide.ns = TRUE, label.y = 40) +
  theme_minimal(base_size = 8) +
  theme(
    axis.text.x = element_text(size = 12, color = "black", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 11, color = "black"),
    axis.title.y = element_text(size = 15, face = "plain"),
    strip.text = element_text(size = 10, face = "bold"),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
    legend.position = "none",
    aspect.ratio = 2,  # ← Keeps it narrow/tall
    plot.margin = margin(4, 4, 4, 4))

p2

# Save as a pdf
pdf("5.1_Myeloid_proportions_3-6mo_allcelltypes.pdf", width = 8, height = 6)
p2
dev.off()


