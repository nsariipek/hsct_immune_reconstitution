# Calculate the myeloid cells proportions
# Load the libraries
library(tidyverse)
library(janitor)
library(ggrepel)
library(ggpubr)

# Empty environment
rm(list=ls())

# Set directory
setwd("~/TP53_ImmuneEscape/5_Cell_Proportions/")

# Load the seurat meta data that I saved previously instead of the whole seurat object
seu_df <- read_csv("~/seu_df_250411.csv")
# seu_df has P32 as relapse sample, so need to change

# Add the TP53 MT information
seu_df <- seu_df %>%
  mutate(TP53_status = case_when(
    patient_id %in% c("P01", "P02", "P03", "P04", "P05", "P06", "P07", "P08", "P09", "P10", "P11", "P12", "P14", "P17") ~ "MT",
    patient_id %in% c("P13", "P15", "P16", "P18", "P19", "P20", "P21", "P22", "P23", "P24", "P25", "P26", "P27",
                      "P28", "P29", "P30", "P31", "P32", "P33") ~ "WT",
    TRUE ~ NA_character_  # fallback in case of unmatched ID
  ))

celltype_levels <- c(
  "Progenitors", "Early Erythroids", "Mid Erythroids", "Late Erythroids",
  "Pro Monocytes", "Monocytes", "Non Classical Monocytes", "cDC", "pDC",
  "Pro B cells", "Pre-B", "B cells", "Plasma cells",
  "CD4 Naïve", "CD4 Effector Memory", "CD4 Memory", "Treg",
  "CD8 Naïve", "CD8 Effector", "CD8 Memory", "CD8 Exhausted",
  "γδ T", "NK T", "Adaptive NK", "CD56 Bright NK", "CD56 Dim NK",
  "Cycling T-NK cells", "UD1", "UD2", "UD3")

# Convert celltype column to a factor with desired order
seu_df$celltype <- factor(seu_df$celltype, levels = celltype_levels)


# Set the sample order same as Figure 1 swimmer plot 
sample_order <- c(paste0("P", str_pad(1:4, 2, pad = "0")),     
                  paste0("P", str_pad(13:27, 2, pad = "0")),     
                  paste0("P", str_pad(5:12, 2, pad = "0")),   
                  paste0("P", str_pad(28:33, 2, pad = "0")))     

bar_data <- seu_df %>%
  mutate(patient_id = factor(patient_id, levels = sample_order)) %>%
  filter(library_type == "MNC",timepoint %in% c("3","5","6") & sample_status =="remission" )%>%
  group_by(patient_id, celltype) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(patient_id) %>%
  mutate(percent = count / sum(count) * 100)

# Load colors from 2.3_PvG-Colors.R
celltype_colors_df <- read.table("../celltype_colors.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE, comment.char = "")
celltype_colors <- setNames(celltype_colors_df$color, celltype_colors_df$celltype)

# Survival colors
survival_colors <- c("Non-relapsed" = "#546fb5FF","Relapsed" = "#e54c35ff")

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
pdf("5.1_Allcells_proportions_3-6mo.pdf", width = 10, height = 5)
p1
dev.off()

################################################################################################################

# Calculate total cells and myeloid cells per patient ==percentage of all myeloid cells out of all cell types in each sample (MNC libraries only)

proportions_df <- seu_df %>%
  filter(timepoint %in% c("3", "5", "6"),
         sample_status == "remission", library_type=="MNC") %>%
  filter(!celltype %in% c("UD1", "UD2","UD3")) %>%
  group_by(patient_id, survival) %>%
  mutate(total_cells = n()) %>%  # total cells per sample
  count(patient_id, survival, total_cells, celltype, name = "cell_count") %>%
  mutate(percent = (cell_count / total_cells) * 100) %>%
mutate(survival = factor(survival, levels = c("Non-relapsed", "Relapsed")))


p2 <- ggplot(proportions_df, aes(x = survival, y = percent, fill = survival)) +
  geom_boxplot(width = 0.7, outlier.shape = NA, alpha = 0.9) +
  geom_jitter(shape = 21, size = 2, stroke = 0.2, color = "black", width = 0.15) +
  facet_wrap(~ celltype, ncol = 9) +
  scale_fill_manual(values = survival_colors) +
  labs(y = "% of all cells", x = NULL) +
  ggpubr::stat_compare_means(
    aes(x = survival, y = percent, group = survival),
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

# Addition if you want to calculate the p-values separately 
library(broom)

pvals_df <- proportions_df %>%
  group_by(celltype) %>%
  filter(!is.na(percent), !is.na(survival)) %>%
  summarise(
    p.value = tryCatch(wilcox.test(percent ~ survival)$p.value, error = function(e) NA_real_),
    .groups = "drop"
  ) %>%
  mutate(label = case_when(
    p.value < 0.001 ~ "***",
    p.value < 0.01  ~ "**",
    p.value < 0.05  ~ "*",
    TRUE ~ "ns"
  ))

