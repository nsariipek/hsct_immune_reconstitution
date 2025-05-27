# Peter van Galen, 250527
# Generate TCR diversity UMAP for Figure 3A, B

# Load libraries
library(tidyverse)
library(Seurat)
library(harmony)
library(RColorBrewer)
library(SeuratWrappers)

# Start with a clean slate
rm(list=ls())

# Set working directory and load Seurat object (Nurefsan, Terra)
#setwd("~/TP53_ImmuneEscape/4_Trajectories/")
#seu <- readRDS("../AuxiliaryFiles/250426_Seurat_annotated.rds")

# Set working directory and load Seurat object (Peter, local)
setwd("~/DropboxMGB/Projects/ImmuneEscapeTP53/TP53_ImmuneEscape/4_Trajectories/")
seu <- readRDS("../AuxiliaryFiles/250426_Seurat_annotated.rds")

# Load colors
celltype_colors_df <- read.table("../celltype_colors.txt", sep = "\t", header = T, stringsAsFactors = F, comment.char = "")
celltype_colors <- setNames(celltype_colors_df$color, celltype_colors_df$celltype)


######## Part 1 - Proccessing using Seurat ########

# Add TCAT annotations
usage_tib <- read_tsv("../3_DGE/3.1_starCAT/starCAT.scores.txt.gz") %>% rename("...1" = "cell")
tcat_celltypes_df <- data.frame(select(usage_tib, Multinomial_Label, Proliferation, Proliferation_binary), row.names = usage_tib$cell)
seu <- AddMetaData(seu, tcat_celltypes_df)
seu$Multinomial_Label <- factor(seu$Multinomial_Label, levels = c("CD4_Naive", "CD4_CM", "CD4_EM", "Treg", "CD8_Naive", "CD8_CM", "CD8_EM",  "CD8_TEMRA", "MAIT", "gdT"))

# Subset for T cells
seu_subset <- subset(seu, Multinomial_Label %in% levels(seu$Multinomial_Label))

# Run standard Seurat steps
seu_subset <- NormalizeData(seu_subset)
seu_subset <- FindVariableFeatures(seu_subset)
seu_subset <- ScaleData(seu_subset)
seu_subset <- RunPCA(seu_subset)

# Run Harmony to remove the batch effect
seu_subset <- RunHarmony(object = seu_subset, group.by.vars = c("patient_id"), plot_convergence = T)
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
seu_subset <- RunUMAP(seu_subset, reduction = "harmony", dims = 1:ndim, return.model = T)

# Visualize the UMAP
p1 <- DimPlot(seu_subset, reduction = "umap", group.by = "Multinomial_Label", label = T, shuffle = T) +
  scale_color_manual(values = celltype_colors) +
  theme(aspect.ratio = 1)
ggsave(paste0("4.2.1_UMAP_Multinomial-Label.pdf"), plot = p1, width = 7, height = 6)


######## Part 2 - Monocle3 Workflow ########

# Convert the Seurat object to the cell_data_set object for monocle3
cds <- as.cell_data_set(seu_subset)

# See cell metadata and counts
#colData(cds)
#counts(cds)[1:20,1:20]

# Assign partitions. This is needed to learn trajectory later
recreate.partition <- as.factor(c(rep(1, length(cds@colData@rownames))))
names(recreate.partition) <- cds@colData@rownames
cds@clusters$UMAP_BMM$partitions <- recreate.partition

# Add cluster info and UMAP coordinates to cell_data_set
cds@clusters$UMAP$clusters <- seu_subset$seurat_clusters
cds@int_colData@listData$reducedDims$UMAP <- seu_subset@reductions$umap@cell.embeddings

# Plot TCAT labels and Seurat clusters
plot_cells(cds, color_cells_by = 'Multinomial_Label',
           label_groups_by_cluster = FALSE,
           group_label_size = 5) +
  theme(legend.position = "right", aspect.ratio = 1)

plot_cells(cds, color_cells_by = 'cluster',
           label_groups_by_cluster = FALSE,
           group_label_size = 5) +
  theme(legend.position = "right", aspect.ratio = 1)


######## Part 3 - Learn trajectory ########

# This takes some time. Do not use partitions because the analysis is already subsetted to closely related T cells
cds <- learn_graph(cds, use_partition = F)

# Select naive T cells as the root, the clusters change depending on the environment
if (T_marker == "CD8") {
  naive_cluster <- 3
} else if (T_marker == "CD4") {
  naive <- c(0,1)
}
cds <- order_cells(cds, reduction_method = "UMAP",
  root_cells = colnames(cds[, clusters(cds) %in% naive]))

plot_cells(cds,
            color_cells_by = "pseudotime",
            label_branch_points = F,
            label_roots = F,
            label_leaves = F) +
  theme(aspect.ratio = 1)
#ggsave...

# Cells ordered by monocle3 pseudotime
cds$monocle3_pseudotime <- pseudotime(cds)
data.pseudo <- as.data.frame(colData(cds))

as.data.frame(colData(cds)) %>%
  mutate(monocle3_pseudotime = pseudotime(cds)) %>%
  ggplot(aes(x = monocle3_pseudotime,
           y = reorder(Multinomial_Label, monocle3_pseudotime, median))) +
  geom_boxplot(outlier.size = 0.8, outlier.alpha = 0.5, width = 0.6) +
  scale_fill_brewer(palette = "Set2") +
  labs(x = "Monocle3 Pseudotime", y = "Cell type") +
  theme_bw() +
  theme(aspect.ratio = 0.5,
    panel.grid = element_blank(),
    legend.position = "none")

ggsave(paste0("4.1_", T_marker, "_pseudotime_boxplot.pdf"), width = 8, height = 6)


######## Part  - Make histograms of pseudotime values by cohorts ########

# Wrangle Seurat object with CD8 T cells
seu_subset$pseudotime <- pseudotime(cds)
Idents(seu_subset) <- seu_subset$celltype

# Some informative visualizations
FeaturePlot(seu_subset, reduction = "umap", features = "pseudotime") + theme(aspect.ratio = 1)
FeaturePlot(seu_subset, reduction = "umap", features = "Proliferation") + theme(aspect.ratio = 1)
DimPlot(seu_subset, reduction = "umap", group.by = "patient_id") + theme(aspect.ratio = 1)

# Extract metadata to facilitate histogram
metadata_tib <- tibble(seu_subset@meta.data, rownames = "cell")

# Subset for time point and mutation status of interest
meta_subset <- metadata_tib %>% filter(sample_status == "remission",
  TP53_status == "MUT", 
  timepoint %in% c(3, 5, 6),
  patient_id != "P21") # Consider excluding P21 T cells - there are very few and they have an outlier cell state distribution

# What is the number of cells per patient?
meta_subset$patient_id %>% table %>% sort %>% rev

# Subset meta_subset for the same number of cells per patient
n_cells <- min(table(meta_subset$patient_id)[table(meta_subset$patient_id) != 0])
meta_subset <- meta_subset %>%
  mutate(patient_id = as.character(patient_id)) %>%
  group_by(patient_id) %>%
  slice_sample(n = n_cells)

#define colors
cohort_colors <- c("long-term-remission" = "#546fb5FF","relapse" = "#e54c35ff")

# Plot histogram per cohort
meta_subset %>%
  ggplot(aes(x = pseudotime, color = cohort)) +
  scale_color_manual(values= cohort_colors)+
  geom_density() +
  theme_bw() +
  theme(aspect.ratio = 0.5, panel.grid = element_blank())

ggsave(paste0("4.1_", T_marker, "_pseudotime_histogram_onlyTP53_MT.pdf"), width = 8, height = 6)

# Plot histogram per patient
meta_subset %>% ggplot(aes(x = pseudotime, color = patient_id)) +
  geom_density() +
  scale_color_manual(values = c(rep("black", 5), rep("red", 4))) +
  theme_bw() +
  theme(aspect.ratio = 0.5, panel.grid = element_blank())

# Statistical test
group1 <- meta_subset %>% filter(cohort == "long-term-remission") %>%
  pull(pseudotime)
group2 <- meta_subset %>% filter(cohort == "relapse") %>%
  pull(pseudotime)
ks.test(group1, group2)

# Violin plot with transparent grey symbols for cells on top
meta_subset %>%
  ggplot(aes(x = cohort, y = pseudotime)) +
  geom_violin() +
  geom_jitter(color = "#80808080")


######## Part 7 - Finding genes that change as a function of pseudotime  ########
 
deg_cd8tcells <- graph_test(cds, neighbor_graph = 'principal_graph', cores = 4)
 
deg_cd8tcells_df <- deg_cd8tcells %>%
  arrange(q_value) %>%
  filter(status == 'OK') %>% head()
 
#write.csv(deg_cd8tcells_df, "deg_cd8tcells_df.csv")
 
FeaturePlot(seu_subset, features = c('TNFRSF18', 'RPL22', 'TNFRSF25'), reduction = "umap")
 
 
 