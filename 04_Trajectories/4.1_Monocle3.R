# Nurefsan Sariipek and Peter van Galen, 250526
# Perform trajectory analysis. The dimensionality reduction and clustering differs slightly depending on the environment - the final version was run by Peter

# Load libraries
library(tidyverse)
library(Seurat)
library(harmony)
library(ggrastr)
library(monocle3)
library(SeuratWrappers)

# Start with a clean slate
rm(list = ls())

# Set working directory and load Seurat object
# fmt: skip
setwd("~/DropboxMGB/Projects/ImmuneEscapeTP53/TP53_ImmuneEscape/4_Trajectories/")
# For VM:
#setwd("~/TP53_ImmuneEscape/4_Trajectories/")
seu <- readRDS("../AuxiliaryFiles/250528_Seurat_complete.rds")


######## Part 1 - Processing using Seurat ########

# Visualize proliferation (filter by <0.05 below)
plot(seu$TCAT_Proliferation, pch = ".")
abline(h = 0.05, col = "red")

# Subset for non-proliferating CD8 T cells
seu_subset <- subset(
  seu,
  subset = TCAT_Multinomial_Label %in%
    c("CD8_Naive", "CD8_CM", "CD8_EM", "CD8_TEMRA") &
    TCAT_Proliferation < 0.05
)
T_marker <- "CD8"
# Or subset for non-proliferating CD4 T cells
seu_subset <- subset(
  seu,
  subset = TCAT_Multinomial_Label %in%
    c("CD4_Naive", "CD4_CM", "CD4_EM") &
    TCAT_Proliferation < 0.05
)
T_marker <- "CD4"

# OPTIONAL: SKIP TO LINE 203

# Run standard Seurat steps
seu_subset <- NormalizeData(seu_subset)
seu_subset <- FindVariableFeatures(seu_subset)
seu_subset <- ScaleData(seu_subset)
seu_subset <- RunPCA(seu_subset)

# Run Harmony to remove the batch effect
seu_subset <- RunHarmony(
  object = seu_subset,
  group.by.vars = c("patient_id"),
  plot_convergence = T
)
ElbowPlot(seu_subset, reduction = "pca")
ElbowPlot(seu_subset, reduction = "harmony")

## Decide on the dimensions by checking different ones
#dims_to_test <- seq(5, 15)

# Create output directory if it doesn't exist
#if (!dir.exists("~/umap_dim_checks")) dir.create("~/umap_dim_checks")

#for (d in dims_to_test) {
#  cat("Running UMAP with dims = 1 :", d, "\n")

#  # Copy the object
#  seu_temp <- seu_subset

#  # Recalculate neighbors and UMAP
#  seu_temp <- FindNeighbors(seu_temp, dims = 1:d, verbose = FALSE)
#  seu_temp <- RunUMAP(seu_temp, reduction = "harmony", dims = 1:d, return.model = TRUE, verbose = FALSE)

#  # Create the plot
#  p <- DimPlot(seu_temp, reduction = "umap", group.by = "TCAT_Multinomial_Label", shuffle = TRUE) +
#    ggtitle(paste0("Dims: 1-", d)) +
#    theme(aspect.ratio = 1)

#  # Save the plot
#  ggsave(
#    filename = paste0("~/umap_dim_checks/dims_", d, ".pdf"),
#    plot = p,
#    width = 6,
#    height = 6
#  )
#}

# Move forward with manually selected number of dimensions - 11 works well for both CD4 and CD8.
ndim <- 11
seu_subset <- FindNeighbors(seu_subset, reduction = "harmony", dims = 1:ndim)
seu_subset <- FindClusters(seu_subset, resolution = 0.5)
seu_subset <- RunUMAP(
  seu_subset,
  reduction = "harmony",
  dims = 1:ndim,
  return.model = T
)

# Visualize the UMAP
p1 <- DimPlot(
  seu_subset,
  reduction = "umap",
  group.by = "TCAT_Multinomial_Label",
  pt.size = 6,
  shuffle = T,
  raster = T,
  raster.dpi = c(1536, 1536)
) +
  theme(aspect.ratio = 1)
DimPlot(
  seu_subset,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = T,
  shuffle = T
) +
  theme(aspect.ratio = 1)
ggsave(
  paste0("4.1.1_", T_marker, "_UMAP_Multinomial-Label.pdf"),
  plot = p1,
  width = 7,
  height = 6
)


######## Part 2 - Monocle3 Workflow ########

# Convert the Seurat object to the cell_data_set object for monocle3
cds <- as.cell_data_set(seu_subset)

# See cell metadata and counts
#colData(cds)
#counts(cds)[1:20,1:20]

# Assign partitions. This is needed to learn trajectory later
recreate.partition <- as.factor(rep(1, length(cds@colData@rownames)))
names(recreate.partition) <- cds@colData@rownames
cds@clusters$UMAP_BMM$partitions <- recreate.partition

# Add cluster info and UMAP coordinates to cell_data_set
cds@clusters$UMAP$clusters <- seu_subset$seurat_clusters
cds@int_colData@listData$reducedDims$UMAP <- seu_subset@reductions$umap@cell.embeddings

# Plot TCAT labels and Seurat clusters
plot_cells(
  cds,
  color_cells_by = 'TCAT_Multinomial_Label',
  label_groups_by_cluster = FALSE,
  group_label_size = 5
) +
  theme(legend.position = "right", aspect.ratio = 1)

plot_cells(
  cds,
  color_cells_by = 'cluster',
  label_groups_by_cluster = FALSE,
  group_label_size = 5
) +
  theme(legend.position = "right", aspect.ratio = 1)

# Learn trajectory: this takes some time. Do not use partitions because the analysis is already subsetted to closely related T cells
cds <- learn_graph(cds, use_partition = F)

# Select naive T cells as the root. The clusters change depending on the environment
if (T_marker == "CD8") {
  naive_clusters <- 3
} else if (T_marker == "CD4") {
  naive_clusters <- c(0, 1)
}
cds <- order_cells(
  cds,
  reduction_method = "UMAP",
  root_cells = colnames(cds[, clusters(cds) %in% naive_clusters])
)

p2 <- plot_cells(
  cds,
  color_cells_by = "pseudotime",
  label_branch_points = F,
  label_roots = F,
  label_leaves = F
) +
  theme(aspect.ratio = 1)
p_raster <- rasterize(p2, dpi = 300)
ggsave(
  paste0("4.1.2_", T_marker, "_Monocle3_pseudotime_UMAP.pdf"),
  plot = p_raster
)


######## Part 3 - Make histograms of pseudotime values by cohorts ########

# Add pseudotime column to Seurat object
seu_subset$monocle3_pseudotime <- pseudotime(cds)

# Alternatively, if skipping Part 2, update seu_subset as follows
saved_data <- read.csv(
  paste0("4.1_", T_marker, "_UMAP_and_pseudotime.csv"),
  row.names = 1
)
seu_subset[["umap"]] <- CreateDimReducObject(
  embeddings = as.matrix(saved_data[, 1:2])
)
seu_subset <- AddMetaData(
  seu_subset,
  metadata = saved_data[, "monocle3_pseudotime", drop = F]
)

# Some informative visualizations
FeaturePlot(seu_subset, reduction = "umap", features = "monocle3_pseudotime") +
  theme(aspect.ratio = 1)
FeaturePlot(seu_subset, reduction = "umap", features = "TCAT_Proliferation") +
  theme(aspect.ratio = 1)
DimPlot(seu_subset, reduction = "umap", group.by = "patient_id") +
  theme(aspect.ratio = 1)

# Extract metadata to facilitate histogram
metadata_tib <- tibble(seu_subset@meta.data, rownames = "cell")

# Subset for time point and mutation status of interest
meta_subset <- metadata_tib %>%
  filter(
    sample_status == "remission",
    timepoint %in% c(3, 5, 6),
    TP53_status == "MUT"
  )

# What is the number of cells per patient?
meta_subset$patient_id %>% table %>% sort %>% rev

# Exclude P21 who has 266 CD4 T cells and 197 CD8 T cells
meta_subset <- meta_subset %>% filter(patient_id != "P21")

# Subset meta_subset for the same number of cells per patient. This leaves us with 913 CD4 T cells or 434 CD8 T cells per patient
n_cells <- min(table(meta_subset$patient_id)[
  table(meta_subset$patient_id) != 0
])
meta_subset <- meta_subset %>%
  mutate(patient_id = as.character(patient_id)) %>%
  group_by(patient_id) %>%
  slice_sample(n = n_cells)

# Boxplot per cell type
x_lim <- c(
  min(meta_subset$monocle3_pseudotime),
  max(meta_subset$monocle3_pseudotime)
)
p3 <- meta_subset %>%
  ggplot(aes(
    x = monocle3_pseudotime,
    y = reorder(TCAT_Multinomial_Label, monocle3_pseudotime, median)
  )) +
  geom_boxplot(outlier.size = 0.8, outlier.alpha = 0.5, width = 0.8) +
  coord_cartesian(xlim = x_lim) +
  labs(y = "Cell type") +
  theme_bw() +
  theme(
    aspect.ratio = 0.2,
    axis.text = element_text(color = "black"),
    axis.ticks = element_line(color = "black"),
    panel.grid = element_blank(),
    legend.position = "none"
  )

# Plot histogram per cohort
p4 <- meta_subset %>%
  ggplot(aes(x = monocle3_pseudotime, color = cohort)) +
  scale_color_manual(
    values = c("long-term-remission" = "#546fb5FF", "relapse" = "#e54c35ff")
  ) +
  coord_cartesian(xlim = x_lim) +
  geom_density() +
  theme_bw() +
  theme(
    aspect.ratio = 0.5,
    axis.text.y = element_text(color = "black"),
    axis.ticks.y = element_line(color = "black"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    panel.grid = element_blank()
  )

p4 / p3 + plot_layout(heights = c(0.8, 0.3))

ggsave(
  paste0("4.1.3_", T_marker, "_boxplot_histogram_TP53-MT_3-6M.pdf"),
  width = 8,
  height = 5
)

# Plot histogram per patient
meta_subset %>%
  ggplot(aes(x = monocle3_pseudotime, color = cohort, linetype = patient_id)) +
  geom_density() +
  scale_color_manual(
    values = c("long-term-remission" = "#546fb5FF", "relapse" = "#e54c35ff")
  ) +
  theme_bw() +
  theme(aspect.ratio = 0.5, panel.grid = element_blank())
ggsave(
  paste0("4.1.4_", T_marker, "_histogram_per-patient.pdf"),
  width = 8,
  height = 6
)

# Statistical test
group1 <- meta_subset %>%
  filter(cohort == "long-term-remission") %>%
  pull(monocle3_pseudotime)
group2 <- meta_subset %>%
  filter(cohort == "relapse") %>%
  pull(monocle3_pseudotime)
ks.test(group1, group2)$p.value

# Save csv file to facilitate reproduction
data2save <- seu_subset@reductions$umap@cell.embeddings
all(rownames(data2save) == colnames(seu_subset))
data2save <- cbind(
  data2save,
  monocle3_pseudotime = seu_subset$monocle3_pseudotime
)
data2save <- apply(data2save, 2, round, 4)
write.csv(
  data2save,
  file = paste0("4.1_", T_marker, "_UMAP_and_pseudotime.csv")
)
