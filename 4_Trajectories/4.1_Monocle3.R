# Nurefsan Sariipek, 250430
# Perform trajectory analysis
# Load the library
library(tidyverse)
library(Seurat)
library(harmony)
library(RColorBrewer)
library(monocle3)
library(SeuratWrappers)

# Start with a clean slate
rm(list=ls())

# Set working directory and load Seurat object (Nurefsan, Terra)
setwd("~/TP53_ImmuneEscape/4_Trajectories/")
seu <- readRDS("~/250428_Tcells.rds")
# # Set working directory and load Seurat object (Peter, local)
# setwd("~/DropboxMGB/Projects/ImmuneEscapeTP53/TP53_ImmuneEscape/2_Annotate-predict/")
# seu <- readRDS("../AuxiliaryFiles/250426_Seurat_annotated.rds")


######## Part 1 - Proccessing using Seurat ########

######## APPROACH 1 (NEW)
    # Add TCAT annotations
    usage_tib <- read_tsv("../3_DGE/3.1_starCAT/starCAT.scores.txt.gz") %>% rename("...1" = "cell")
    tcat_celltypes_df <- data.frame(select(usage_tib, Multinomial_Label, Proliferation, Proliferation_binary), row.names = usage_tib$cell)
    seu <- AddMetaData(seu, tcat_celltypes_df)
    seu$Multinomial_Label <- factor(seu$Multinomial_Label, levels = c("CD4_Naive", "CD4_CM", "CD4_EM", "Treg", "CD8_Naive", "CD8_CM", "CD8_EM",  "CD8_TEMRA", "MAIT", "gdT"))
    plot(usage_tib$Proliferation, pch = ".")
    abline(h = 0.05, col = "red")

    # Subset for CD8 T cells
    seu_subset <- subset(seu, subset = Multinomial_Label %in% c("CD8_Naive", "CD8_CM", "CD8_EM",  "CD8_TEMRA") & Proliferation < 0.05)

    # Or subset for CD4 T cells
    seu_subset <- subset(seu, subset = Multinomial_Label %in% c("CD4_Naive", "CD4_CM", "CD4_EM") & Proliferation < 0.05)

######## APPROACH 2 (OLD)

    # Load the seurat object that contains the T cells only

    # Subset for CD8 T cells
    seu_subset <- subset(seu, subset = celltype %in% c("CD8 Naive","CD8 Central Memory","CD8 Effector Memory 1","CD8 Effector Memory 2"))
    #"CD8 Tissue Resident Memory", "T Proliferating"

    # Subset for CD4 T cells (should also change cd8 to cd4 below)
    seu_subset <- subset(seu, subset = celltype %in% c("CD4 Central Memory", "CD4 Naive" ,"CD4 Effector Memory","CD4 Regulatory" ))

# Normalize
seu_subset <- NormalizeData(seu_subset)
seu_subset <- FindVariableFeatures(seu_subset)

# Scale the data (commented out cell cycle regression b/c Proliferation_binary is excluded)
#seu_subset <- CellCycleScoring(seu_subset,
#  s.features = cc.genes.updated.2019$s.genes,
#  g2m.features = cc.genes.updated.2019$g2m.genes)
#seu_subset <- ScaleData(seu_subset, vars.to.regress = c("S.Score", "G2M.Score"))
seu_subset <- ScaleData(seu_subset)
seu_subset <- RunPCA(seu_subset)

# Run Harmony to remove the batch effect
ElbowPlot(seu_subset)
seu_subset <- RunHarmony(object = seu_subset, group.by.vars = c("patient_id"), plot_convergence = T)
ElbowPlot(seu_subset, reduction = "harmony")

## Decide on the dimensions by checking different ones
dims_to_test <- seq(5, 15)

# Create output directory if it doesn't exist
if (!dir.exists("~/umap_dim_checks")) dir.create("~/umap_dim_checks")

for (d in dims_to_test) {
  cat("Running UMAP with dims = 1 :", d, "\n")

  # Copy the object
  seu_temp <- seu_subset

  # Recalculate neighbors and UMAP
  seu_temp <- FindNeighbors(seu_temp, dims = 1:d, verbose = FALSE)
  seu_temp <- RunUMAP(seu_temp, reduction = "harmony", dims = 1:d, return.model = TRUE, verbose = FALSE)

  # Create the plot
  p <- DimPlot(seu_temp, reduction = "umap", group.by = "Multinomial_Label", shuffle = TRUE) +
    ggtitle(paste0("Dims: 1-", d)) +
    theme(aspect.ratio = 1)

  # Save the plot immediately
  ggsave(
    filename = paste0("~/umap_dim_checks/dims_", d, ".pdf"),
    plot = p,
    width = 6,
    height = 6
  )
}

# Move forward with manually selected number of dimensions, for CD4 12 looks the best
seu_subset <- FindNeighbors(seu_subset, reduction = "harmony", dims = 1:12)
seu_subset <- FindClusters(seu_subset, resolution = 1)
seu_subset <- RunUMAP(seu_subset, reduction = "harmony", dims = 1:12, return.model = T)

# Visualize the UMAP
DimPlot(seu_subset, reduction = "umap", group.by = "Multinomial_Label", shuffle = T) +
  theme(aspect.ratio = 1)
DimPlot(seu_subset, reduction = "umap", group.by = "seurat_clusters", label = T, shuffle = T) +
  theme(aspect.ratio = 1)


######## Part 2 - Monocle3 Workflow ########

# Convert the Seurat object to the cell_data_set object for monocle3
cds <- as.cell_data_set(seu_subset)

# See cell metadata and gene metadata
#colData(cds)
#fData(cds)
#rownames(fData(cds))[1:10]

# Add gene_short_name column
#fData(cds)$gene_short_name <- rownames(fData(cds))

# See counts
#counts(cds)[1:20,1:20]

######## Part 3 - Add clusters from Seurat ########

# Assign partitions
recreate.partition <- c(rep(1, length(cds@colData@rownames)))
names(recreate.partition) <- cds@colData@rownames
recreate.partition <- as.factor(recreate.partition)

cds@clusters$UMAP_BMM$partitions <- recreate.partition

# Assign the cluster info
cds@clusters$UMAP$clusters <- seu_subset$seurat_clusters

# Assign UMAP coordinate - cell embeddings
cds@int_colData@listData$reducedDims$UMAP <- seu_subset@reductions$umap@cell.embeddings

# Plot clusters
plot_cells(cds, color_cells_by = 'cluster',
           label_groups_by_cluster = FALSE,
           group_label_size = 5) +
  theme(legend.position = "right", aspect.ratio = 1)


plot_cells(cds, color_cells_by = 'Multinomial_Label',
           label_groups_by_cluster = FALSE,
           group_label_size = 5) +
  #scale_color_manual(values= c("red","blue","green","maroon","yellow","gray","cyan"))+
  theme(legend.position = "right", aspect.ratio = 1)



######## Part 4 - Learn trajectory ########

# Do not use partitions because the analysis is already subsetted to closely related CD8+ T cells
cds_temp <- learn_graph(cds, use_partition = F)

plot_cells(cds_temp,
           color_cells_by = "celltype",
           label_branch_points = FALSE,
           label_roots = FALSE,
           label_leaves = FALSE,
           group_label_size = 5)

######## Part 5 - Order the cells in pseudotime ########

# Select naive T cells as the root, the clusters change for Peter and Nurefsan, this is for Nurefsan
# CD4 0,4,6,9,11
# CD8 3,13
cds_temp <- order_cells(cds_temp, reduction_method = "UMAP",
  root_cells = colnames(cds_temp[, clusters(cds_temp) %in% c(0,4,6,9,11)]))

plot_cells(cds_temp,
            color_cells_by = "pseudotime",
            label_groups_by_cluster = FALSE,
            label_branch_points = FALSE,
            label_roots = FALSE,
            label_leaves = FALSE,
            group_label_size = 5)

# Cells ordered by monocle3 pseudotime
cds_temp$monocle3_pseudotime <- pseudotime(cds_temp)
data.pseudo <- as.data.frame(colData(cds_temp))
 
ggplot(data.pseudo, 
       aes(x = monocle3_pseudotime, 
           y = reorder(Multinomial_Label, monocle3_pseudotime, median), 
           fill = Multinomial_Label)) +
  geom_boxplot(outlier.size = 0.8, outlier.alpha = 0.5, width = 0.6) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    x = "Monocle3 Pseudotime",
    y = "Cell Type (ordered by median pseudotime)",
    fill = "Cell Type"
  ) +
  theme_classic(base_size = 13) +  # Classic theme = no grid lines
  theme(
    axis.text.x = element_text(size = 11, color = "black"),
    axis.text.y = element_text(size = 11, color = "black"),
    axis.title = element_text(size = 13),
    axis.line = element_line(size = 0.8),
    legend.position = "none"
  )

ggsave("4.1_CD8_pseudotime_boxplot.pdf", width = 8, height = 6)
ggsave("4.1_CD4_pseudotime_boxplot.pdf", width = 8, height = 6)
######## Part 6 - Make histograms of pseudotime values by cohorts ########

# Wrangle Seurat object with CD8 T cells
seu_subset$pseudotime <- pseudotime(cds_temp)
Idents(seu_subset) <- seu_subset$celltype

# Some informative visualizations
FeaturePlot(seu_subset, reduction = "umap", features = "pseudotime", split.by = "cohort") + theme(aspect.ratio = 1)
ggsave("4.1_CD8_pseudotime_umap_cohort.pdf", width = 8, height = 6)
ggsave("4.1_CD4_pseudotime_umap_cohort.pdf", width = 8, height = 6)
FeaturePlot(seu_subset, reduction = "umap", features = "Proliferation") + theme(aspect.ratio = 1)
DimPlot(seu_subset, reduction = "umap", group.by = "celltype") + theme(aspect.ratio = 1)
DimPlot(seu_subset, reduction = "umap", group.by = "Multinomial_Label") + theme(aspect.ratio = 1)
DimPlot(seu_subset, reduction = "umap", group.by = "patient_id") + theme(aspect.ratio = 1)

# Extract metadata to facilitate histogram
metadata_tib <- tibble(seu_subset@meta.data, rownames = "cell")
metadata_tib$umap_1 <- seu_subset@reductions$umap@cell.embeddings[,1]
metadata_tib$umap_2 <- seu_subset@reductions$umap@cell.embeddings[,2]

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

ggsave("4.1_CD4_pseudotime_histogram_onlyTP53_MT.pdf", width = 8, height = 6)
ggsave("4.1_CD8_pseudotime_histogram_onlyTP53_MT.pdf", width = 8, height = 6)

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
 
deg_cd8tcells <- graph_test(cds_temp, neighbor_graph = 'principal_graph', cores = 4)
 
deg_cd8tcells_df <- deg_cd8tcells %>%
  arrange(q_value) %>%
  filter(status == 'OK') %>% head()
 
#write.csv(deg_cd8tcells_df, "deg_cd8tcells_df.csv")
 
FeaturePlot(seu_subset, features = c('TNFRSF18', 'RPL22', 'TNFRSF25'), reduction = "umap")
 
 
 