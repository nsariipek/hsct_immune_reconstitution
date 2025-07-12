# Nurefsan Sariipek and Peter van Galen, 250712
# Visualize the Numbat results

# Load Libraries
library(tidyverse)
library(Seurat)
library(janitor)
#library(readr)
#library(RColorBrewer)
#library(ggsci)
#library(scattermore)
#library(ggh4x)

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
  header = TRUE,
  stringsAsFactors = FALSE,
  comment.char = ""
)
celltype_colors <- setNames(
  celltype_colors_df$color,
  celltype_colors_df$celltype
)
compartment_colors <- c(tumor = "#4D4D4D", normal = "#D9B88C")
timepoint_colors <- c(`pre-transplant` = "#A3BFD9", relapse = "#8B0000")

#### Plot Numbat compartment ####
DimPlot(seu_numbat, group.by = "numbat_compartment", shuffle = T) +
  scale_color_manual(values = compartment_colors) +
  theme(
    aspect.ratio = 1,
    axis.line = element_blank(),
    panel.border = element_rect(color = "black")
  )

# Save
ggsave("8.5.1_UMAP_genotype.pdf", width = 7, height = 7)

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
ggsave("8.5.2_UMAP_status.pdf", width = 15, height = 7)

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
ggsave("8.5.3_UMAP_celltype.pdf", width = 15, height = 7)

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

ggsave("8.5.4_Numbat_compartment_per_status.pdf", width = 8, height = 4)

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

ggsave("8.5.5_Celltype_barplots.pdf", width = 12, height = 6)

#### Pseudotime of tumor cells ####

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
meta_subset <- meta_subset %>%
  group_by(patient_id, sample_status) %>%
  slice_sample(n = 127)

meta_subset %>%
  ggplot(aes(x = predicted_Pseudotime, color = sample_status)) +
  geom_density(linewidth = 2) +
  scale_color_manual(values = timepoint_colors) +
  theme_bw() +
  theme(aspect.ratio = 0.3, panel.grid = element_blank())

ggsave("8.5.6_Pseudotime_maligant_cells.pdf", width = 8, height = 6)

# Statistical test
group1 <- meta_subset %>%
  filter(sample_status == "pre-transplant") %>%
  pull(predicted_Pseudotime)
group2 <- meta_subset %>%
  filter(sample_status == "relapse") %>%
  pull(predicted_Pseudotime)
ks.test(group1, group2)

# Boxplot per cell type (to place below the pseudotime plot)
meta_subset %>%
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

ggsave("8.5.7_Pseudotime_celltypes.pdf", width = 12, height = 4)