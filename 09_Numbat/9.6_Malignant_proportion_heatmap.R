# Peter van Galen, 250705
# Plot the proportion of malignant cells per cell type (x) and patient (y) as a heatmap

# Load libraries
library(Seurat)
library(tidyverse)
library(janitor)

# Set working directory
repo_root <- system("git rev-parse --show-toplevel", intern = TRUE)
setwd(paste0(repo_root, "/09_Numbat"))

# Clear environment variables
rm(list = ls())

# Load data
seu <- readRDS("../AuxiliaryFiles/250528_Seurat_complete.rds")

# Filter, then count number of cells per patient_id, celltype, and compartment
cell_counts <- as_tibble(seu@meta.data) %>%
  filter(
    !is.na(numbat_compartment),
    !is.na(celltype)
  ) %>%
  dplyr::count(patient_id, celltype, numbat_compartment)

# Pivot wider to compute proportion malignant
malignant_props <- cell_counts %>%
  pivot_wider(
    names_from = numbat_compartment,
    values_from = n,
    values_fill = 0
  ) %>%
  mutate(proportion_malignant = tumor / (tumor + normal))

# Plot
heatmap <- malignant_props %>%
  mutate(patient_id = factor(patient_id, levels = rev(unique(patient_id)))) %>%
  complete(patient_id, celltype) %>%
  ggplot(aes(x = celltype, y = patient_id, fill = proportion_malignant)) +
  geom_raster() +
  scale_fill_gradientn(
    colors = c(
      "#2166AC",
      "#4393C3",
      "#92C5DE",
      "#F4A582",
      "#D6604D",
      "#B2182B"
    ),
    limits = c(0, 1),
    name = "Percentage",
    na.value = "grey80"
  ) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(
    y = "Patient ID",
    title = "Proportion malignant cells (of recipient cells, all time points)"
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, ),
    axis.title.x = element_blank(),
    axis.ticks = element_line(color = "black")
  ) +
  coord_fixed(ratio = 1)

# View
heatmap

pdf("9.6_Malignant_proportion_heatmap.pdf", width = 9, height = 4)
heatmap
dev.off()


# For the legend:
as_tibble(seu@meta.data) %>%
  filter(
    !is.na(numbat_compartment),
    !is.na(celltype)
  ) %>%
  dplyr::count(patient_id, sample_status) %>%
  pivot_wider(
    names_from = sample_status,
    values_from = n
  )
# --> "Cells from all time points are included and P33 mainly represents pre-transplant cells."
