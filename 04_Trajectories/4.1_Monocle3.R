# Nurefsan Sariipek and Peter van Galen, 250526
# Perform trajectory analysis. The dimensionality reduction and clustering differs slightly depending on the environment - the final version was run by Peter

# Load libraries
library(tidyverse)
library(Seurat)
library(SeuratWrappers)
library(harmony)
library(ggrastr)
library(monocle3)
library(patchwork)
library(ggpubr)

# Start with a clean slate
rm(list = ls())

######## LOAD DATA AND SUBSET TO RELEVANT CELL TYPE ########

# Set working directory and load Seurat object
# fmt: skip
setwd("~/DropboxMGB/Projects/ImmuneEscapeTP53/hsct_immune_reconstitution/04_Trajectories/")
# For VM:
#setwd("~/hsct_immune_reconstitution/4_Trajectories/")
seu <- readRDS("../AuxiliaryFiles/250528_Seurat_complete.rds")

# Load colors
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

# Visualize proliferation (filter by <0.05 below)
plot(seu$TCAT_Proliferation, pch = ".")
abline(h = 0.05, col = "red")

# Subset for non-proliferating CD4 T cells
seu_subset <- subset(
  seu,
  subset = TCAT_Multinomial_Label %in%
    c("CD4_Naive", "CD4_CM", "CD4_EM") &
    TCAT_Proliferation < 0.05
)
T_marker <- "CD4"
# Or subset for non-proliferating CD8 T cells
seu_subset <- subset(
  seu,
  subset = TCAT_Multinomial_Label %in%
    c("CD8_Naive", "CD8_CM", "CD8_EM", "CD8_TEMRA") &
    TCAT_Proliferation < 0.05
)
T_marker <- "CD8"

# OPTIONAL: SKIP TO LINE 215

######## GENERATE UMAPS USING SEURAT ########

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

######## MONOCLE3 WORKFLOW ########

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

######## PSEUDOTIME DENSITY PLOT ########

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

# Determine bin breaks (this becomes more relevant in the stats section)
bin_breaks <- seq(
  from = min(meta_subset$monocle3_pseudotime),
  to = max(meta_subset$monocle3_pseudotime),
  length.out = 4
)

# Boxplot per cell type
x_lim <- c(
  min(meta_subset$monocle3_pseudotime),
  max(meta_subset$monocle3_pseudotime)
)
p3 <- meta_subset %>%
  ggplot(aes(
    x = monocle3_pseudotime,
    y = reorder(TCAT_Multinomial_Label, monocle3_pseudotime, median),
    fill = TCAT_Multinomial_Label
  )) +
  geom_boxplot(outlier.size = 0.8, outlier.alpha = 0.5, width = 0.8) +
  scale_fill_manual(values = celltype_colors) +
  coord_cartesian(xlim = x_lim) +
  labs(y = "TCAT Label") +
  theme_bw() +
  theme(
    aspect.ratio = 0.1,
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
  geom_vline(
    xintercept = bin_breaks,
    linetype = "dashed",
    alpha = 0.7,
    color = "black"
  ) +
  theme_bw() +
  theme(
    aspect.ratio = 0.25,
    axis.text.y = element_text(color = "black"),
    axis.ticks.y = element_line(color = "black"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    panel.grid = element_blank()
  )

p4 / p3 + plot_layout(heights = c(0.2, 0.1))

ggsave(
  paste0("4.1.3_", T_marker, "_boxplot_histogram_TP53-MT_3-6M.pdf"),
  width = 8,
  height = 3
)

# Plot histogram per patient
meta_subset %>%
  ggplot(aes(x = monocle3_pseudotime, color = cohort, linetype = patient_id)) +
  geom_density() +
  geom_vline(
    xintercept = bin_breaks,
    alpha = 0.7,
    color = "black"
  ) +
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

######## STATISTICS ########

# We initially thought about comparing pseudotime values between the cohorts using the Kolmogorov-Smirnov test or Wilcoxon test. However, using cells as replicates is not OK.
# Instead, we should compare patient medians. In this section, we will compare the proportion of cells that fall within three pseudotime bins.

# Create 3 equal-width bins of pseudotime values
meta_binned <- meta_subset %>%
  select(patient_id, cohort, TCAT_Multinomial_Label, monocle3_pseudotime) %>%
  mutate(
    pseudotime_bin = cut(
      monocle3_pseudotime,
      breaks = 3,
      labels = c("Low", "Medium", "High")
    )
  )

# Calculate proportions of cells per patient/cohort in each bin
proportions_data <- meta_binned %>%
  group_by(patient_id, cohort, pseudotime_bin) %>%
  summarize(count = n()) %>%
  mutate(proportion = count / sum(count)) %>%
  ungroup()

proportions_data %>%
  ggplot(aes(x = cohort, y = proportion, fill = cohort)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter() +
  scale_fill_manual(
    values = c("long-term-remission" = "#546fb5FF", "relapse" = "#e54c35ff")
  ) +
  facet_wrap(~pseudotime_bin) +
  stat_compare_means(
    method = "wilcox.test",
    label = "p.format",
    label.x.npc = 0.5,
    label.y.npc = 0.95,
    hjust = 0.5
  ) +
  theme_bw() +
  theme(
    aspect.ratio = 2,
    axis.text = element_text(color = "black"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank()
  )

ggsave(
  paste0("4.1.5_", T_marker, "_bin_proportions_per-patient.pdf"),
  width = 7,
  height = 4
)


######## SAVE ########

# Save cell embeddings and pseudotime values  to facilitate reproduction
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
