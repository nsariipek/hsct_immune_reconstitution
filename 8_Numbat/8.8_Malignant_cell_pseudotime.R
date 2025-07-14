# Nurefsan Sariipek and Peter van Galen, 250712
# Visualize the Numbat results

# Load libraries
library(tidyverse)
library(Seurat)
library(janitor)

# Empty environment
rm(list = ls())

# Set working directory (for Nurefsan)
setwd("~/TP53_ImmuneEscape/8_Numbat/")

# For Peter
setwd("~/DropboxMGB/Projects/ImmuneEscapeTP53/TP53_ImmuneEscape/8_Numbat")

# Load the saved Seurat object
seu <- readRDS("../AuxiliaryFiles/250528_Seurat_complete.rds")

# Subset for cells with Numbat calls
seu_numbat <- subset(seu, !is.na(numbat_compartment))

# Define colors to use in plots
celltype_colors_df <- read.table(
  "../celltype_colors.txt",
  sep = "\t",
  header = T,
  stringsAsFactors = F,
  comment.char = ""
)
celltype_colors <- setNames(
  celltype_colors_df$color,
  celltype_colors_df$celltype
)
timepoint_colors <- c(`pre-transplant` = "#A3BFD9", relapse = "#8B0000")

# Extract metadata
metadata_tib <- tibble(seu_numbat@meta.data, rownames = "cell")

# Subset for malignant cell types of interest
metadata_tib %>%
  filter(numbat_compartment == "tumor") %>%
  dplyr::count(celltype) %>%
  print(n = 40)

meta_subset <- metadata_tib %>%
  filter(
    numbat_compartment == "tumor",
    celltype %in%
      c(
        "HSC MPP",
        "MEP",
        "EoBasoMast Precursor",
        "Megakaryocyte Precursor",
        "Early Erythroid",
        "Late Erythroid",
        "LMPP",
        "Cycling Progenitor",
        "Early GMP",
        "Late GMP",
        "Pro-Monocyte",
        "CD14 Mono",
        "CD16 Mono",
        "cDC",
        "pDC"
      )
  )

# What samples can we use?
meta_subset %>%
  dplyr::count(patient_id, sample_status) %>%
  pivot_wider(names_from = sample_status, values_from = n) %>%
  adorn_totals(where = "row")
# There are not enough remission cells to include this analysis, and we should only use samples with pre-tx and relapse time points
meta_subset <- meta_subset %>%
  filter(sample_status != "remission", !patient_id %in% c("P22", "P33"))

# Subset for the same number of cells per patient
meta_subset %>% dplyr::count(patient_id, sample_status)
meta_subset2 <- meta_subset %>%
  group_by(patient_id, sample_status) %>%
  slice_sample(n = 127)

meta_subset2 %>%
  ggplot(aes(x = predicted_Pseudotime, color = sample_status)) +
  geom_density(linewidth = 2) +
  scale_color_manual(values = timepoint_colors) +
  theme_bw() +
  theme(aspect.ratio = 0.3, panel.grid = element_blank())

ggsave("8.8.1_Pseudotime_maligant_cells.pdf", width = 8, height = 4)

# Statistical test
group1 <- meta_subset2 %>%
  filter(sample_status == "pre-transplant") %>%
  pull(predicted_Pseudotime)
group2 <- meta_subset2 %>%
  filter(sample_status == "relapse") %>%
  pull(predicted_Pseudotime)
ks.test(group1, group2)

# Boxplot per cell type (to place below the pseudotime plot)
meta_subset2 %>%
  ggplot(aes(x = predicted_Pseudotime, y = 1, fill = celltype)) +
  geom_tile(width = 0.1, linewidth = 0) +
  scale_fill_manual(values = celltype_colors) +
  labs(y = "Cell type") +
  theme_bw() +
  theme(
    aspect.ratio = 0.05,
    axis.text.x = element_text(color = "black"),
    axis.ticks.x = element_line(color = "black"),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    legend.position = "bottom"
  )

ggsave("8.8.2_Pseudotime_celltypes.pdf", width = 12, height = 4)

