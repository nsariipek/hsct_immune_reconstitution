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
compartment_colors <- c(normal = "#D9B88C", tumor = "#4D4D4D")
timepoint_colors <- c(`pre-transplant` = "#A3BFD9", relapse = "#8B0000")


#### Plot Numbat compartment by sample status ####
DimPlot(
  seu_numbat,
  group.by = "numbat_compartment",
  shuffle = T,
  split.by = "sample_status",
  pt.size = 1,
  alpha = 0.8
) +
  scale_color_manual(values = compartment_colors) +
  theme(
    aspect.ratio = 1,
    axis.line = element_blank(),
    panel.border = element_rect(color = "black")
  )

# Annotate in Illustrator
seu_numbat@meta.data %>% dplyr::count(sample_status, numbat_compartment)

# Save
ggsave("8.7.1_UMAP_status.pdf", width = 15, height = 7)


#### Plot celltype by Numbat compartment ####
DimPlot(
  subset(seu_numbat, !is.na(celltype)),
  group.by = "celltype",
  shuffle = T,
  split.by = "numbat_compartment"
) +
  scale_color_manual(values = celltype_colors) +
  theme(
    aspect.ratio = 1,
    axis.line = element_blank(),
    panel.border = element_rect(color = "black")
  )

# Save
ggsave("8.7.2_UMAP_celltype.pdf", width = 15, height = 7)


#### Normal vs. maligant composition per status ####
seu_numbat@meta.data %>%
  dplyr::count(patient_id, numbat_compartment, sample_status) %>%
  ggplot(aes(x = patient_id, y = n, fill = numbat_compartment)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = compartment_colors) +
  facet_wrap(~sample_status) +
  ylab("Fraction of cells") +
  xlab("Patient") +
  ggtitle("Normal vs. maligant composition per status") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave("8.7.3_Numbat_compartment_per_status.pdf", width = 8, height = 4)

#### Plot cell type composition per patient ####
df_prop <- seu_numbat@meta.data %>%
  dplyr::count(cohort, patient_id, numbat_compartment, celltype) %>%
  group_by(cohort, patient_id, numbat_compartment) %>%
  mutate(prop = n / sum(n))

ggplot(df_prop, aes(x = numbat_compartment, y = prop, fill = celltype)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = celltype_colors) +
  ylab("Proportion of cells") +
  xlab("") +
  facet_wrap(~patient_id, nrow = 1) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave("8.7.4_Celltype_barplots.pdf", width = 12, height = 6)
