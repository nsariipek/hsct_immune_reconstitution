# Nurefsan Sariipek and Peter van Galen, 250715
# Calculate CD4/8 Ratio

# Load libraries
library(tidyverse)
library(Seurat)
library(ggpubr)

# Set working directory (for Nurefsan):
setwd("~/hsct_immune_reconstitution/5_Cell_Proportions/")

# For Peter
# fmt: skip
setwd("~/DropboxMGB/Projects/ImmuneEscapeTP53/hsct_immune_reconstitution/5_Cell_Proportions")

# Clear environment variables
rm(list = ls())

# Load Seurat object & subset for T cells
seu <- readRDS("../AuxiliaryFiles/250528_Seurat_complete.rds")
seu_T <- subset(seu, !is.na(TCAT_Multinomial_Label))

# Cohort colors
cohort_colors <- c("long-term-remission" = "#546fb5", "relapse" = "#e54c35")

# Make a dataframe and add CD4/CD8
metadata_df <- seu_T@meta.data %>%
  dplyr::select(
    patient_id,
    cohort,
    timepoint,
    sample_status,
    TP53_status,
    celltype
  )
metadata_df$type <- case_when(
  grepl("CD8.", metadata_df$celltype) ~ "CD8",
  grepl("CD4.", metadata_df$celltype) ~ "CD4"
)
metadata_df$type <- factor(metadata_df$type, levels = c("CD4", "CD8"))

# Select only 100-day samples that were in remission at that timepoint
meta_subset <- metadata_df %>%
  filter(
    sample_status == "remission",
    timepoint %in% c(3, 5, 6),
    TP53_status == "MUT" # Change to see WT samples or remove to see all samples
  )

# Calculate the ratio for each patient
tb <- meta_subset %>%
  group_by(patient_id, cohort, type) %>%
  count() %>%
  ungroup() %>%
  group_by(patient_id) %>%
  pivot_wider(names_from = type, values_from = n) %>%
  mutate(ratio = CD4 / CD8)

# Look at results
tb
tb %>% group_by(cohort) %>% summarize(median(ratio))

# Plot
p1 <- tb %>%
  ggplot(aes(x = cohort, y = ratio, fill = cohort)) +
  geom_boxplot(width = 0.7, alpha = 0.9, outlier.shape = NA) +
  geom_jitter(shape = 21, size = 2, color = "black") +
  scale_fill_manual(values = cohort_colors) +
  coord_cartesian(ylim = c(0, 2)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  labs(y = "CD4/CD8 T cell ratio", x = NULL) +
  stat_compare_means(
    aes(group = cohort),
    method = "wilcox.test",
    label = "p.format",
    label.y = 2,
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

# View
p1

# Save as pdf (depending on filter)
pdf("5.3_CD4-CD8_ratio_All-boxplot.pdf", width = 2, height = 3)
p1
dev.off()

pdf("5.3_CD4-CD8_ratio_MUT-boxplot.pdf", width = 2, height = 3)
p1
dev.off()

pdf("5.3_CD4-CD8_ratio_WT-boxplot.pdf", width = 2, height = 3)
p1
dev.off()
