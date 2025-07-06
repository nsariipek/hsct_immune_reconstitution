# Peter van Galen, 250705
# Plot the proportion of malignant cells per cell type (x) and patient (y) as a heatmap

# Load libraries
library(Seurat)
library(tidyverse)
library(janitor)

# Set working directory
setwd("~/DropboxMGB/Projects/ImmuneEscapeTP53/TP53_ImmuneEscape/8_Numbat")

# Delete environment variables
rm(list = ls())

# Load data
seu <- readRDS("../AuxiliaryFiles/250528_Seurat_complete.rds")

# Filter, then count number of cells per patient_id, celltype, and compartment
cell_counts <- as_tibble(seu@meta.data) %>%
  filter(
    !is.na(numbat_compartment),
    !is.na(celltype),
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
  labs(
    x = "Cell Type",
    y = "Patient ID",
    title = "Malignant cell proportion (all time points)"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_text(
      size = 8,
      angle = 45,
      hjust = 1,
      vjust = 1,
      color = "black"
    ),
    axis.text.y = element_text(size = 8, color = "black"),
    axis.title.x = element_text(size = 10, color = "black"),
    axis.title.y = element_text(size = 11, color = "black"),
    plot.title = element_text(size = 12, color = "black", hjust = 0.5),
    legend.title = element_text(size = 11, color = "black"),
    panel.grid = element_blank(),
    axis.ticks = element_line(color = "black", size = 0.5)
  ) +
  coord_fixed(ratio = 1)

pdf("8.7_Malignant_proportion_heatmap.pdf", width = 8, height = 8)
heatmap
dev.off()
