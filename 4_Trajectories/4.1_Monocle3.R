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
setwd("~/TP53_ImmuneEscape/2_Annotate-predict/")
seu <- readRDS("~/250428_Tcells.rds")
# Set working directory and load Seurat object (Peter, local)
setwd("~/DropboxMGB/Projects/ImmuneEscapeTP53/TP53_ImmuneEscape/2_Annotate-predict/")
seu <- readRDS("../AuxiliaryFiles/250426_Seurat_annotated.rds")


######## Part 1 - Proccessing using Seurat ########

######## APPROACH 1 (NEW)
    # Add TCAT annotations
    usage_tib <- read_tsv("../3_DGE/3.1_starCAT/starCAT.scores.txt.gz") %>% rename("...1" = "cell")
    tcat_celltypes_df <- data.frame(select(usage_tib, Multinomial_Label), row.names = usage_tib$cell)
    seu <- AddMetaData(seu, tcat_celltypes_df)
    seu$Multinomial_Label <- factor(seu$Multinomial_Label, levels = c("CD4_Naive", "CD4_CM", "CD4_EM", "Treg", "CD8_Naive", "CD8_CM", "CD8_EM",  "CD8_TEMRA", "MAIT", "gdT"))

    # Subset for CD8 T cells
    seu_cd8 <- subset(seu, subset = Multinomial_Label %in% c("CD8_Naive", "CD8_CM", "CD8_EM",  "CD8_TEMRA"))


######## APPROACH 2 (OLD)

    # Load the seurat object that contains the T cells only

    # Subset for CD8 T cells
    seu_cd8 <- subset(seu, subset = celltype %in% c("CD8 Naive","CD8 Central Memory","CD8 Effector Memory 1","CD8 Effector Memory 2"))
    #"CD8 Tissue Resident Memory", "T Proliferating"

    # Subset for CD4 T cells (should also change cd8 to cd4 below)
    seu_cd4 <- subset(seu, subset = celltype %in% c("CD4 Central Memory", "CD4 Naive" ,"CD4 Effector Memory","CD4 Regulatory" ))

# Normalize
seu_cd8 <- NormalizeData(seu_cd8)
seu_cd8 <- FindVariableFeatures(seu_cd8)
# Scale the data
seu_cd8 <- ScaleData(seu_cd8)
seu_cd8 <- RunPCA(seu_cd8)
# Run Harmony to remove the batch effect
ElbowPlot(seu_cd8)
seu_cd8 <- RunHarmony(object = seu_cd8, group.by.vars = c("patient_id"), plot_convergence = T)
ElbowPlot(seu_cd8, reduction = "harmony")

# ## Decide on the dimensions by checking different ones
# dims_to_test <- seq(10, 20, by = 2)
# 
# # Create output directory if it doesn't exist
# if (!dir.exists("umap_dim_checks")) dir.create("umap_dim_checks")
# 
# for (d in dims_to_test) {
#   cat("Running UMAP with dims = 1:", d, "\n")
# 
#   # Copy the object
#   seu_temp <- seu_cd8
# 
#   # Recalculate neighbors and UMAP
#   seu_temp <- FindNeighbors(seu_temp, dims = 1:d, verbose = FALSE)
#   seu_temp <- RunUMAP(seu_temp, reduction = "harmony", dims = 1:d, return.model = TRUE, verbose = FALSE)
# 
#   # Create the plot
#   p <- DimPlot(seu_temp, reduction = "umap", group.by = "celltype", shuffle = TRUE) +
#     ggtitle(paste("Dims: 1-", d)) +
#     theme(aspect.ratio = 1)
# 
#   # Save the plot immediately
#   ggsave(
#     filename = paste0("umap_dim_checks/dims_", d, ".pdf"),
#     plot = p,
#     width = 6,
#     height = 6
#   )
# }

# Move forward with manually selected number of dimensions
seu_cd8 <- FindNeighbors(seu_cd8, reduction = "harmony", dims = 1:10)
seu_cd8 <- RunUMAP(seu_cd8, reduction = "harmony", dims = 1:10, return.model = T)

# Visualize the UMAP
DimPlot(seu_cd8, reduction = "umap", group.by = "Multinomial_Label", shuffle = T) +
  theme(aspect.ratio = 1)


######## Part 2 - Monocle3 Workflow ########

# Convert the seurat object to the cell_data_set object for monocle3
cds <- as.cell_data_set(seu_cd8)
# See cell metadata
colData(cds)
# See gene metadata
fData(cds)
rownames(fData(cds))[1:10]
# Add gene_short_name column
fData(cds)$gene_short_name <- rownames(fData(cds))
# See counts
counts(cds)[1:20,1:20]

######## Part 3 - Cluster cells using clustering info from Seurat's UMAP ########

# Assign partitions
reacreate.partition <- c(rep(1,length(cds@colData@rownames)))
names(reacreate.partition) <- cds@colData@rownames
reacreate.partition <- as.factor(reacreate.partition)

cds@clusters$UMAP_BMM$partitions <- reacreate.partition

# Assign the cluster info
#list_cluster <- seu_cd8@active.ident
list_cluster <- seu_cd8$Multinomial_Label
cds@clusters$UMAP$clusters <- list_cluster

# Assign UMAP coordinate - cell embeddings
cds@int_colData@listData$reducedDims$UMAP <- seu_cd8@reductions$umap@cell.embeddings

# Plot clusters
cluster.before.trajectory <- plot_cells(cds,
           color_cells_by = 'cluster',
           label_groups_by_cluster = FALSE,
           group_label_size = 5)+
  theme(legend.position = "right")

cluster.before.trajectory

cluster.names <- plot_cells(cds, color_cells_by = 'Multinomial_Label',
                                       label_groups_by_cluster = FALSE,
                                       group_label_size = 5) +
  #scale_color_manual(values= c("red","blue","green","maroon","yellow","gray","cyan"))+
  theme(legend.position = "right")

cluster.names


######## Part 4 - Learn trajectory ########

# Not completely sure if use_partition = T is the right choice. This takes some time.
cds_temp <- learn_graph(cds, use_partition = TRUE)

plot_cells(cds_temp,
           color_cells_by = "celltype",
           label_branch_points = FALSE,
           label_roots = FALSE,
           label_leaves = FALSE,
           group_label_size = 5)

######## Part 5 - Order the cells in pseudotime ########

# Select CD8 Naive cells in the last part (clusters 9 and 11 for Nurefsan, 8 and 10 for Peter)
cds_temp <- order_cells(cds_temp, reduction_method = "UMAP", root_cells = colnames(cds_temp[, clusters(cds_temp) %in% c("CD8_Naive")]))

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
  aes(monocle3_pseudotime, reorder(celltype, monocle3_pseudotime, median),
  fill= celltype)) +
  geom_boxplot()


######## Part 6 - Make histograms of pseudotime values by cohorts ########

# Wrangle Seurat object with CD8 T cells
seu_cd8$pseudotime <- pseudotime(cds_temp)
Idents(seu_cd8) <- seu_cd8$celltype

# Some informative visualizations
FeaturePlot(seu_cd8, reduction= "umap", features = "pseudotime", label = T, split.by = "cohort")
DimPlot(seu_cd8, reduction= "umap", group.by = "celltype")
DimPlot(seu_cd8, reduction= "umap", group.by = "patient_id")
DimPlot(seu_cd8, reduction= "umap", group.by = "Multinomial_Label")
 
# Extract metadata to facilitate histogram
metadata_tib <- tibble(seu_cd8@meta.data, rownames = "cell")
metadata_tib$umap_1 <- seu_cd8@reductions$umap@cell.embeddings[,1]
metadata_tib$umap_2 <- seu_cd8@reductions$umap@cell.embeddings[,2]

# Subset for time point and mutation status of interest
meta_subset <- metadata_tib %>% filter(sample_status == "remission",
  TP53_status == "MUT", timepoint %in% c(3, 5, 6),
  patient_id != "P21") # TEST/TMP

# What is the number of cells per patient?
meta_subset$patient_id %>% table %>% sort %>% rev

# Subset meta_subset for the same number of cells per patient
n_cells <- min(table(meta_subset$patient_id)[table(meta_subset$patient_id) != 0])
meta_subset <- meta_subset %>%
  mutate(patient_id = as.character(patient_id)) %>%
  group_by(patient_id) %>%
  slice_sample(n = n_cells)

# Plot histogram per cohort
meta_subset %>% ggplot(aes(x = pseudotime, color = cohort)) +
  geom_density(bw = 0.5) +
  theme_bw() +
  theme(aspect.ratio = 0.5, panel.grid = element_blank())

# Plot histogram per patient
meta_subset %>% ggplot(aes(x = pseudotime, color = patient_id)) +
  geom_density() +
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
 
FeaturePlot(seu_cd8, features = c('TNFRSF18', 'RPL22', 'TNFRSF25'), reduction = "umap")
 
 
 