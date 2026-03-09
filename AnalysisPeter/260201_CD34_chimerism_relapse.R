# Peter van Galen, 260201
# Plot UMAP colored by CD34 expression and outcome correlation

library(tidyverse)
library(Seurat)
library(ggrastr)
library(ggrepel)
library(patchwork)

# Set working directory
repo_root <- system("git rev-parse --show-toplevel", intern = TRUE)
setwd(paste0(repo_root, "/AnalysisPeter"))

# Clear environment variables & load favorite function
rm(list = ls())
cutf <- function(x, f = 1, d = "/") {
  sapply(strsplit(x, d), function(i) paste(i[f], collapse = d))
}

# Load data
seu <- readRDS("../AuxiliaryFiles/250528_Seurat_complete.rds")
seu <- NormalizeData(seu)

# Extract metadata from Seurat object
metadata_tib <- as_tibble(seu@meta.data)

# Add some columns
metadata_tib$UMAP_1 <- seu@reductions$umap_bmm@cell.embeddings[, 1]
metadata_tib$UMAP_2 <- seu@reductions$umap_bmm@cell.embeddings[, 2]
metadata_tib$CD34 <- LayerData(seu, layer = "data")["CD34", ]


# Relate chimerism to cell type -----------------------------------------------

# Subset data to relevant time points (similar to 08_Souporcell/8.4_Souporcell_plots.R)
metasubset_tib <- metadata_tib |>
  filter(
    !is.na(celltype),
    souporcell_origin %in% c("donor", "recipient"),
    sample_status == "remission",
    timepoint %in% c(3, 5, 6)
  )

# Calculate donor chimerism for each cohort/cell type, and the difference between cohorts. Note that this could be dominated by patients with more cells since we're summarizing by cohort
celltype_chimerism_diff <- metasubset_tib |>
  group_by(celltype, cohort) |>
  summarize(donor_prop = sum(souporcell_origin == "donor") / n()) |>
  pivot_wider(names_from = cohort, values_from = donor_prop) |>
  mutate(prop_diff = `long-term-remission` - relapse)

# Join data for plotting
plot_tib <- left_join(
  select(metasubset_tib, celltype, UMAP_1, UMAP_2, CD34),
  select(celltype_chimerism_diff, celltype, prop_diff)
)

# Visualize
p1 <- plot_tib |>
  mutate(mean_CD34 = mean(CD34), .by = "celltype") |>
  ggplot(aes(x = UMAP_1, y = UMAP_2, color = mean_CD34)) +
  geom_point_rast(size = 0.3) +
  scale_color_gradientn(
    colors = c("grey90", "grey90", "red"),
    values = scales::rescale(c(0, 0.5, 1.5))
  ) +
  theme_bw() +
  theme(
    aspect.ratio = 1,
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  )

p2 <- plot_tib |>
  ggplot(aes(x = UMAP_1, y = UMAP_2, color = prop_diff)) +
  geom_point_rast(size = 0.3) +
  scale_color_gradientn(
    colors = c("grey90", "grey90", "red"),
    values = scales::rescale(c(0, 0.3, 0.8))
  ) +
  theme_bw() +
  theme(
    aspect.ratio = 1,
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  )

# Plot
p1 + p2

# Correlation
cor_data <- plot_tib |>
  mutate(mean_CD34 = mean(CD34), .by = "celltype") |>
  select(celltype, mean_CD34, prop_diff) |>
  unique()

p3 <- cor_data %>%
  ggplot(aes(x = mean_CD34, y = prop_diff)) +
  geom_point(color = "red", size = 2) +
  geom_text_repel(
    aes(label = ifelse(mean_CD34 > 0.35, as.character(celltype), "")),
    size = 3,
    max.overlaps = 20
  ) +
  labs(
    x = "Mean CD34 expression",
    y = "Chimerism difference (LTR - Relapse)"
  ) +
  theme_bw() +
  theme(
    aspect.ratio = 1,
    axis.text = element_text(color = "black"),
    axis.ticks = element_line(color = "black"),
    panel.grid = element_blank(),
  )

# Test (to be added in Illustrator)
cor.test(cor_data$mean_CD34, cor_data$prop_diff)

p1 + p2 + p3

# Save
ggsave("260201_CD34_chimerism_relapse.pdf", width = 12, height = 4.5)
