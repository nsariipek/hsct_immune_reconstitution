# Peter van Galen, 250527
# Generate TCR diversity UMAP for Figure 3A, B

# Load libraries
library(tidyverse)
library(Seurat)
library(harmony)
library(RColorBrewer)
library(SeuratWrappers)
library(viridis)

# Start with a clean slate
rm(list = ls())

# Set working directory (Nurefsan, Terra)
setwd("~/TP53_ImmuneEscape/4_Trajectories/")

# For Peter
# fmt: skip
setwd("~/DropboxMGB/Projects/ImmuneEscapeTP53/TP53_ImmuneEscape/4_Trajectories/")

# Load data
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

# Subset for T cells
seu_subset <- subset(
  seu,
  TCAT_Multinomial_Label %in% levels(seu$TCAT_Multinomial_Label)
)

# TO SKIP SEURAT PROCESSING (WHICH VARIES PER ENVIRONMENT), SKIP TO LINE 109

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
#dims_to_test <- seq(10, 20)

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
#  p <- DimPlot(seu_temp, reduction = "umap", group.by = "Multinomial_Label", shuffle = TRUE) +
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

# Move forward with manually selected number of dimensions. This is optimized for Peter's environment
ndim <- 17
seu_subset <- FindNeighbors(seu_subset, reduction = "harmony", dims = 1:ndim)
#seu_subset <- FindClusters(seu_subset, resolution = 0.5)
seu_subset <- RunUMAP(
  seu_subset,
  reduction = "harmony",
  dims = 1:ndim,
  return.model = T
)

# Save UMAP coordinates to csv file
write.csv(
  as.data.frame(seu_subset@reductions$umap@cell.embeddings),
  file = "4.2_UMAP-embeddings.csv"
)

# Add UMAP coordinates from saved csv file
umap_embeddings <- read.csv(file = "4.2_UMAP-embeddings.csv") %>%
  column_to_rownames(var = "X")
seu_subset[["umap"]] <- CreateDimReducObject(as.matrix(umap_embeddings))

# Visualize the UMAP
p1 <- DimPlot(
  seu_subset,
  reduction = "umap",
  group.by = "TCAT_Multinomial_Label",
  raster = T,
  raster.dpi = c(1536, 1536),
  pt.size = 3,
  label = T,
  shuffle = T
) +
  scale_color_manual(values = celltype_colors) +
  theme(aspect.ratio = 1)

# View
p1

# Save
ggsave(
  paste0("4.2.1_UMAP_Multinomial-Label.pdf"),
  plot = p1,
  width = 7,
  height = 6
)

# Subset for T cells with TCR info
cells_with_TCR <- rownames(seu_subset@meta.data[!is.na(seu$CTstrict), ])
seu_TCR <- subset(seu_subset, cells = cells_with_TCR)

# Determine clone size
metadata_tib <- as_tibble(seu_TCR@meta.data, rownames = "cell")
metadata_tib <- metadata_tib %>%
  group_by(patient_id, CTstrict) %>%
  reframe(cell = cell, n = n())
metadata_df <- data.frame(metadata_tib, row.names = "cell")
seu_TCR <- AddMetaData(seu_TCR, select(metadata_df, n))

# Cap at 100
seu_TCR$n[seu_TCR$n > 100] <- 100

# Plot TCR diversity
p2 <- FeaturePlot(
  seu_TCR,
  reduction = "umap",
  features = "n",
  raster = T,
  raster.dpi = c(1536, 1536),
  pt.size = 3
) +
  scale_color_viridis_c() +
  theme(aspect.ratio = 1)

# View
p2

# Save
ggsave(paste0("4.2.2_CloneSize.pdf"), plot = p2, width = 7, height = 6)

# Plot TCR diversity for TP53 mutated patients in remission at 3-6M after transplant split by cohort
seu_TCR_subset <- subset(
  seu_TCR,
  sample_status == "remission" &
    timepoint %in% c(3, 5, 6) &
    TP53_status == "MUT"
)

rem <- FeaturePlot(
  subset(seu_TCR_subset, cohort == "long-term-remission"),
  reduction = "umap",
  features = "n",
  raster = T,
  raster.dpi = c(1536, 1536),
  pt.size = 5
) +
  scale_color_viridis_c() +
  theme(aspect.ratio = 1)
rel <- FeaturePlot(
  subset(seu_TCR_subset, cohort == "relapse"),
  reduction = "umap",
  features = "n",
  raster = T,
  raster.dpi = c(1536, 1536),
  pt.size = 5
) +
  scale_color_viridis_c() +
  theme(aspect.ratio = 1)

# View
rem + rel

ggsave(
  paste0("4.2.3_CloneSize_TP53-3-6M_split.pdf"),
  plot = rem + rel,
  width = 7,
  height = 10
)
