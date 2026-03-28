# Peter van Galen, 260328
# Plot pseudotime vs. expression for gene of interest

library(tidyverse)
library(Seurat)
library(ggsci)
library(ggnewscale)

# Set working directory
repo_root <- system("git rev-parse --show-toplevel", intern = TRUE)
setwd(paste0(repo_root, "/11_Expression_across_datasets"))

# Clear environment variables & load favorite function
rm(list = ls())

# Load Seurat data from three studies. This is saved in Terra's Google bucket and available upon request (see 11.1_Merge_and_annotate_three_studies.R).
seu <- readRDS("../AuxiliaryFiles/260110_Seurat_ThreeStudies.rds")

# Combine healthy/malignant calls from different projects
seu$combined_call <- seu@meta.data |>
  mutate(
    combined_call = case_when(
      cell_origin == "Healthy" ~ "healthy", # aml_2019 and healthy_donors_2023
      PredictionRefined == "malignant" ~ "malignant", # aml_2019
      numbat_compartment == "tumor" ~ "malignant" # hsct_2026
    )
  )

# Subset to HSPC/myeloid/erythroid lineages
seu_subset <- subset(seu, celltype %in% c(
  "HSC MPP", "MEP", "EoBasoMast Precursor", "Megakaryocyte Precursor",
  "Early Erythroid", "Late Erythroid", "LMPP", "Cycling Progenitor",
  "Early GMP", "Late GMP", "Pro-Monocyte", "CD14 Mono",
  "CD16 Mono", "cDC", "pDC"
))

# Select gene of interest and add expression values
mygene <- "EIF4A1"
seu_subset$expr <- GetAssayData(seu_subset, layer = "data")[mygene, ]

# Plot all cells
seu_subset@meta.data |>
  filter(combined_call %in% c("healthy", "malignant")) |>
  ggplot(aes(x = predicted_Pseudotime, y = expr, color = combined_call)) +
  geom_point(alpha = 0.1, size = 0.5) +
  geom_smooth() +
  scale_color_manual(
    values = c("healthy" = "#4682B4", "malignant" = "#C0392B")
  ) +
  ylab(paste(mygene, "expression")) +
  theme_bw() +
  theme(aspect.ratio = 2/3,
    panel.grid = element_blank())

ggsave(paste0("11.3_All_cells_", mygene, ".pdf"), width = 8, height = 5)

# Average per donor per cell type
donor_celltype_expr_tib <- seu_subset@meta.data |>
  filter(combined_call %in% c("healthy", "malignant")) |>
  group_by(donor_id, celltype, combined_call) |>
  filter(n() >= 10) |>
  summarize(
    expr = mean(expr),
    pseudotime = mean(predicted_Pseudotime),
    .groups = "drop"
  )

# Statistical test
donor_celltype_expr_tib |>
  summarize(
    p_value = wilcox.test(expr[combined_call == "healthy"], expr[combined_call == "malignant"])$p.value,
    n_healthy = sum(combined_call == "healthy"),
    n_malignant = sum(combined_call == "malignant"),
    .by = celltype
  ) |>
  mutate(p_adj = p.adjust(p_value, method = "BH")) |>
  arrange(p_adj)

# Plot
donor_celltype_expr_tib |>
  ggplot(aes(x = pseudotime, y = expr)) +
  geom_smooth(aes(color = combined_call)) +
  scale_color_manual(
    values = c("healthy" = "#4682B4", "malignant" = "#C0392B")
  ) +
  new_scale_color() +
  geom_point(
    aes(color = donor_id),
    size = 2,
    alpha = 0.7
  ) +
  scale_color_igv(guide = guide_legend(ncol = 2)) +
  ylab(paste(mygene, "expression")) +
  theme_bw() +
  theme(aspect.ratio = 2/3,
    panel.grid = element_blank())

ggsave(paste0("11.3_Averages_", mygene, ".pdf"), width = 8, height = 5)