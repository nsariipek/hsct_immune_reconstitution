# Peter van Galen, 250718
# Evaluate clonotype dominance in remission samples ~3 months after transplant

# Load the libraries
library(tidyverse)
library(Seurat)
library(ggrastr)
library(ggpubr)

# Empty environment
rm(list = ls())

# Set working directory
# fmt: skip
setwd("~/DropboxMGB/Projects/ImmuneEscapeTP53/hsct_immune_reconstitution/06_TCR_Diversity")

# Load final Seurat object including TCR calls & subset for T cells
seu <- readRDS("../AuxiliaryFiles/250528_Seurat_complete.rds")
seu_T <- subset(seu, subset = !is.na(TCAT_Multinomial_Label))

# Cohort colors
cohort_colors <- c("long-term-remission" = "#546fb5", "relapse" = "#e54c35")

# Turn metadata to a tibble and keep only needed variables
metadata_tib <- as_tibble(seu_T@meta.data, rownames = "cell") %>%
  select(
    cell,
    cohort,
    patient_id,
    timepoint,
    sample_status,
    TP53_status,
    CTstrict
  )

# Filter data for remission samples ~3M after transplant
metasubset1_tib <- metadata_tib %>%
  filter(
    sample_status == "remission",
    timepoint %in% c(3, 5, 6)
  )

# With our withoug filter for TP53-mutated patients: choose one!
metasubset2_tib <- metasubset1_tib
patient_filter <- "both"
# Or:
metasubset2_tib <- filter(metasubset1_tib, TP53_status == "MUT")
patient_filter <- "MT"

# View proportion that each clone represents of its sample
metasubset2_tib %>%
  group_by(patient_id, cohort, CTstrict) %>%
  count() %>%
  group_by(patient_id, cohort) %>%
  mutate(proportion = n / sum(n)) %>%
  arrange(desc(proportion)) %>%
  group_by(cohort) %>%
  mutate(sorted_clonotypes = row_number()) %>%
  ggplot(aes(x = sorted_clonotypes, y = proportion, color = cohort)) +
  geom_point_rast(size = 0.5) +
  scale_color_manual(values = cohort_colors) +
  scale_x_log10() +
  scale_y_log10(labels = scales::label_number()) +
  theme_bw() +
  theme(
    aspect.ratio = 1,
    panel.grid = element_blank(),
    axis.text = element_text(color = "black"),
    panel.border = element_rect(color = "black", fill = NA),
    axis.ticks = element_line(color = "black")
  )

# Save
ggsave(
  paste0("6.3.1_Clonotype_proportion_", patient_filter, ".pdf"),
  width = 4.5,
  height = 3
)

# Evaluate the number of clonotypes per patient after downsampling to 300 cells

# Cell threshold
min_cells <- 300

# Subset to samples with a reasonable number of cells
metasubset2_tib %>%
  group_by(patient_id) %>%
  count() %>%
  ggplot(aes(x = patient_id, y = n)) +
  geom_point() +
  geom_hline(yintercept = min_cells, color = "red") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )

# Filter for samples with sufficient cells
metasubset3_tib <- metasubset2_tib %>%
  group_by(patient_id) %>%
  filter(n() >= min_cells)

# Sample min_cells from each group
set.seed(94)
metasubset4_tib <- metasubset3_tib %>%
  group_by(patient_id) %>%
  slice_sample(n = min_cells)

# Plot number of clonotypes in the subset of T cells
metasubset4_tib %>%
  group_by(patient_id, cohort) %>%
  summarize(clonotype_number = length(unique(CTstrict))) %>%
  ggplot(aes(x = cohort, y = clonotype_number, fill = cohort)) +
  geom_boxplot(width = 0.7, alpha = 0.9, outlier.shape = NA) +
  geom_jitter(width = 0.3, shape = 21, size = 2, color = "black") +
  coord_cartesian(ylim = c(0, 310)) +
  scale_fill_manual(values = cohort_colors) +
  stat_compare_means(
    aes(group = cohort),
    method = "wilcox.test",
    label = "p.format",
    label.y = 300
  ) +
  theme_bw() +
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

# Note: the difference for patients with TP53-mutated AML depends on set.seed/slice_sample above

# Save
ggsave(
  paste0("6.3.2_Clonotype_number_", patient_filter, ".pdf"),
  width = 3,
  height = 3.5
)
