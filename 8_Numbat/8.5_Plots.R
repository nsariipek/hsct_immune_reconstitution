# Nurefsan Sariipek, 241219
# Combine Numbat seurat object for UMAP visulazation
# Visualize the Numbat Results, 250401, NS

# Load Libraries
library(readr)
library(Seurat)
library(tidyverse)
library(RColorBrewer)
library(ggsci)
library(scattermore)
library(ggh4x)

# Empty environment
rm(list=ls())

# Set working directory
setwd("~/TP53_ImmuneEscape/8_Numbat/")

# Load the seurat object that has Numbat results from 8.4
seu_combined <- readRDS("~/250505_numbat_combined_seurat.rds")

# Set sample status order
seu_combined$sample_status <- factor(
  seu_combined$sample_status,
  levels = c("pre-transplant", "remission", "relapse"))

seu_combined$compartment_opt <- factor(seu_combined$compartment_opt,levels = c("normal", "tumor"))


# Remove the early relapse patients
seu_combined_subset <- subset(seu_combined, subset=cohort_detail=="1-Relapse")

# Define colors to use in plots
celltype_colors_df <- read.table("../celltype_colors.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE, comment.char = "")
celltype_colors <- setNames(celltype_colors_df$color, celltype_colors_df$celltype)
compartment_colors <- c("tumor" = "#4D4D4D", "normal" = "#FFECB3")
survival_strip_colors <-c("long-term-remission" = "#546fb5FF","relapse" = "#e54c35ff")

# Plot 1
#Visualize the UMAP
UMAP_genotype <- seu_combined_subset@meta.data %>%
  sample_frac(1) %>%
  ggplot(aes(x = UMAP_1, y = UMAP_2, color = compartment_opt)) +
  geom_scattermore(pointsize = 8, pixels = c(4096, 4096)) +
  scale_color_manual(values = compartment_colors)+
  theme_bw() +
  theme(
    aspect.ratio = 1,
    legend.title = element_blank(),
    legend.position = "right") +
  ggtitle("Tumor vs Normal Cells")  

# Save
pdf("UMAP_genotype.pdf", width = 7, height = 7)  
UMAP_genotype
dev.off()

# Plot 2


UMAP_status <- seu_combined_subset@meta.data %>%
  drop_na(UMAP_1, UMAP_2, compartment_opt, sample_status) %>%
  sample_frac(1) %>%
  ggplot(aes(x = UMAP_1, y = UMAP_2, color = compartment_opt)) +
  geom_scattermore(pointsize = 6, pixels = c(2048, 2048)) +
  scale_color_manual(values = compartment_colors) +
  facet_wrap(~sample_status, nrow = 1) +
  theme_bw(base_size = 14) +
  theme(
    aspect.ratio = 1,
    strip.text = element_text(size = 14),
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    strip.background = element_blank(), 
    legend.title = element_blank(),
    legend.position = "right") +
  ggtitle("Tumor vs Normal Cells Across Sample Status")

# Save
pdf("8.5_Relapse_cohort_only_UMAP_status.pdf", width = 7, height = 7)  
UMAP_status
dev.off()

# Plot 3
UMAP_celltype <- seu_combined_subset@meta.data %>%
  drop_na(UMAP_1, UMAP_2, celltype, compartment_opt) %>%
  sample_frac(1) %>%
  ggplot(aes(x = UMAP_1, y = UMAP_2, color = celltype)) +
  geom_scattermore(pointsize = 6, pixels = c(2048, 2048)) +
  scale_color_manual(values = celltype_colors) +
  facet_wrap(~compartment_opt, nrow = 1) +
  theme_bw(base_size = 14) +
  theme(aspect.ratio = 1,
    strip.text = element_text(size = 14),
    strip.background = element_blank(),
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.title = element_blank(),
    legend.position = "none"
  ) 
pdf("UMAP_celltype.pdf", width = 7, height = 7)  
UMAP_celltype
dev.off()

# Plot 4
p4 <- seu_combined_subset@meta.data %>%
  count(patient_id, compartment_opt, sample_status) %>%  # include sample_status
  ggplot(aes(x = patient_id, y = n, fill = compartment_opt)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = compartment_colors) +
  facet_wrap(~sample_status) +
  ylab("Fraction of Cells") +
  xlab("Patient") +
  theme_bw(base_size = 14) +
  ggtitle("Tumor vs Normal Composition per Patient")

pdf("8.5_Relapse_cohort_Numbat_celltype.pdf", width = 8, height = 4)  
p4
dev.off()



# Plot 5
df_prop <- seu_combined_subset@meta.data %>%
  count(cohort, patient_id, compartment_opt, celltype) %>%
  group_by(cohort, patient_id, compartment_opt) %>%
  mutate(prop = n / sum(n))

# Plot
p5 <- ggplot(df_prop, aes(x = compartment_opt, y = prop, fill = celltype)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_grid2(
    rows = vars(cohort),
    cols = vars(patient_id),
    switch = "y",
    strip = strip_themed(background_y = elem_list_rect(fill = survival_strip_colors),
      text_y = elem_list_text(face = "bold", size = 12, color = "white")
    )
  ) +
  scale_fill_manual(values = celltype_colors) +
  ylab("Proportion of Cells (Normalized to 1)") +
  xlab("") +
  theme_bw(base_size = 14) +
  theme(
    strip.placement = "outside",
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_blank(),
    panel.grid = element_blank(),
    panel.spacing = unit(1, "lines")
  )

pdf("8.5_Celltype_Barplots.pdf", width = 16, height = 12)  
p5
dev.off()

