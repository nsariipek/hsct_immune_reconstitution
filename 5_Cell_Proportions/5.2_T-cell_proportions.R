# Nurefsan Sariipek and Peter van Galen, 250715
# Visualize T cell proportions

# Load the libraries
library(tidyverse)
library(Seurat)
library(ggpubr)

# Set working directory (for Nurefsan):
setwd("~/TP53_ImmuneEscape/5_Cell_Proportions/")

# For Peter
# fmt: skip
setwd("~/DropboxMGB/Projects/ImmuneEscapeTP53/TP53_ImmuneEscape/5_Cell_Proportions")

# Clear environment variables
rm(list = ls())

# Load the Seurat object
seu <- readRDS("../AuxiliaryFiles/250528_Seurat_complete.rds")

# Celltype colors
celltype_colors_df <- read.table(
  "../celltype_colors.txt",
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE,
  comment.char = ""
)
celltype_colors <- setNames(
  celltype_colors_df$color,
  celltype_colors_df$celltype
)

# Cohort colors
cohort_colors <- c("long-term-remission" = "#546fb5FF", "relapse" = "#e54c35ff")

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

# Prepare the data
metadata_df <- seu_T@meta.data
celltype_proportions <- metadata_df %>%
  mutate(patient_id = factor(patient_id)) %>%
  filter(
    sample_status == "remission",
    timepoint %in% c(3, 5, 6)
  ) %>%
  group_by(patient_id, celltype) %>%
  summarize(count = n(), .groups = "drop") %>%
  group_by(patient_id) %>%
  mutate(percent = count / sum(count) * 100)

# Plot
p1 <- celltype_proportions %>%
  ggplot(aes(x = patient_id, y = percent, fill = celltype)) +
  geom_bar(stat = "identity", width = 0.8) +
  labs(x = "Patient ID", y = "Percent of T cells") +
  scale_fill_manual(values = celltype_colors) +
  theme_bw() +
  theme(
    aspect.ratio = 0.6,
    panel.grid = element_blank(),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    axis.line = element_blank(),
    axis.ticks = element_line(color = "black")
  )

# View
p1

# Save as a pdf
pdf("5.2_T-cell_proportions_3-6mo_bars.pdf", width = 10, height = 5)
p1
dev.off()

#### Boxplots

# Prepare the data
celltype_proportions2 <- metadata_df %>%
  mutate(patient_id = factor(patient_id)) %>%
  filter(
    timepoint %in% c(3, 5, 6),
    sample_status == "remission",
    TP53_status == "MUT" # Change to see WT samples or remove to see all samples
  ) %>%
  group_by(patient_id, celltype, cohort) %>%
  summarize(count = n(), .groups = "drop") %>%
  group_by(patient_id) %>%
  mutate(percent = count / sum(count) * 100)

# Box plots with P values
p2 <- celltype_proportions2 %>%
  ggplot(aes(x = cohort, y = percent, fill = cohort)) +
  geom_boxplot(width = 0.7, alpha = 0.9, outlier.shape = NA) +
  geom_jitter(shape = 21, size = 2, color = "black") +
  facet_wrap(~celltype, ncol = 10) +
  scale_fill_manual(values = cohort_colors) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  labs(y = "Percent of T cells", x = NULL) +
  stat_compare_means(
    aes(x = cohort, y = percent, group = cohort),
    method = "wilcox.test",
    label = "p.format",
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

# Check the plot
p2

# Save as pdf (depending on filter)
pdf("5.2_T-cell_proportions_All-boxplots.pdf", width = 8, height = 6)
p2
dev.off()

pdf("5.2_T-cell_proportions_TP53_MUT-boxplots.pdf", width = 8, height = 6)
p2
dev.off()

pdf("5.2_T-cell_proportions_TP53_WT-boxplots.pdf", width = 8, height = 6)
p2
dev.off()
