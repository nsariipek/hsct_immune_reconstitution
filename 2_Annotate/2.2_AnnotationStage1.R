# Nurefsan Sariipek and Peter van Galen, updated 250410
# This script quantifies signatures from previous papers (saved in Signatures) and visualizes them on UMAPs (saved in FeaturePlots). In the end, Seurat cluster annotations for the cells (subset 20%) and saved as a csv file

library(tidyverse)
library(Seurat)
library(randomcoloR)
library(R.utils)
#devtools::install_github('immunogenomics/presto')

# Start with a clean slate
rm(list=ls())

setwd("~/TP53_ImmuneEscape/2_Annotate/")

# Load the Seurat object from 2.1_Subset20percent_and_cluster.R that contains 20% of the cells
seu20 <- readRDS("~/250410_SubsettedSeuratObject.rds")


# VISUALIZE MARKER GENES AND SIGNATURES ----------------------------------------

# Visualize the current clusters and save as PDF
p1 <- DimPlot(seu20, reduction = "umap", group.by = "seurat_clusters", shuffle = T, label = T, raster = T) +
  theme(aspect.ratio = 1, legend.position = "none")

pdf("2.2.1_UMAP_clusters.pdf", width = 10, height = 8)
print(p1)
dev.off()

# Find and save the most differentially expressed genes
seu_markers <- FindAllMarkers(seu20, min.pct = .3, logfc.threshold = .3)
top50_markers_tib <- as_tibble(seu_markers) %>% filter(avg_log2FC > 0) %>%
  group_by(cluster) %>% arrange(-avg_log2FC) %>% slice_head(n = 50) %>%
  select(cluster, gene) %>%
  pivot_wider(names_from = "cluster", values_from = "gene",
              names_prefix = "cluster_", values_fn = list) %>%
  unnest(cols = everything())
write_csv(top50_markers_tib, file = "2.2_Stage1_MarkerGenes.csv")

# Plot features to get an idea of cluster identities
FeaturePlot(seu20, features = c("CD34","MPO", "CD14", "MS4A1", "CD3G","CD8B"))
# T cell features
FeaturePlot(seu20, features = c("CD8A", "TCF7", "TOX", "HAVCR2", "CXCR3",
                              "SLAMF6", "CD3E", "CD4", "SELL", "CD44", "PDCD1",
                              "FOXP3", "GZMB", "GZMK", "LAG3", "CD101",
                              "CXCR5", "KLRG1", "IFNG", "TNF")) 
# HSPC markers
FeaturePlot(seu20, features = c("VIM", "FLT3", "CD34", "ITGAL", "THY1",
                              "PTPRC", "KIT", "SLAMF1", "MME", "SLAMF2", "MPO")) 
# B cell markers
FeaturePlot(seu20, features = c("CD19", "MS4A1", "PDCD1LG2", "NT5E", "FCER2", "SDC1",
                              "PAX5", "TCF3", "CD80", "SPIB", "BCL6")) 

# We used previously annotated single cell analysis from our lab for manual annotation, including unpublished data from a discovery cohort:
discoverygenes <- read.table(file = "~/TP53_ImmuneEscape/2_Annotate/Signatures/DiscoveryCohort_MarkerGenes.txt", 
                             sep = "\t",header = TRUE, quote = "", stringsAsFactors = FALSE)
# Add module scores
for (n in names(discoverygenes)) {
  print(n)
  seu20 <- AddModuleScore(object = seu20, features = discoverygenes[n], name = n)
  colnames(seu20@meta.data) <- gsub(str_c(n, "1$"), str_c(n, "_DC_Score"), colnames(seu20@meta.data))
}

# Function to split features into chunks
split_features <- function(features, chunk_size) {
  split(features, ceiling(seq_along(features) / chunk_size))
}

# Define your specific feature groups
feature_groups <- list(
  Group1 = c("CD4.Memory_DC_Score", "CD4.Naïve_DC_Score", "CD56.Bright.NK.cells_DC_Score", 
             "CD56.Dim.NK.cells_DC_Score", "Treg_DC_Score", "γδ.T.lymphocytes_DC_Score", 
             "CD8.Central.Memory_DC_Score", "CD8.Effector.Memory_DC_Score", "CD8.Naïve_DC_Score", 
             "NK.T.cells_DC_Score", "CD8.Terminally.Exhausted_DC_Score"),
  Group2 = c("Doublets_DC_Score", "Early.Erythroids_DC_Score", "pDC_DC_Score", 
             "Late.Erythroids_DC_Score", "Mid.Erythroids_DC_Score", "cDC_DC_Score"),
  Group3 = c("Monocytes_DC_Score", "Non.Classical.Monocytes_DC_Score", 
             "Pro.Monocytes_DC_Score", "Undetermined_DC_Score", "HSPCs_DC_Score"),
  Group4 = c("Pre.B.cells_DC_Score", "Pro.B.cells_DC_Score", 
             "B.cells_DC_Score", "Plasma.Cells_DC_Score"))

# Set the number of features per page
features_per_page <- 4

# Loop through each group and save plots
for (group_name in names(feature_groups)) {
  group_features <- feature_groups[[group_name]]
  
  # Split features into smaller chunks
  feature_chunks <- split_features(group_features, features_per_page)
  
  # Open a PDF file for the group
  pdf(file = paste0("FeaturePlots/", group_name, "_Discovery_FeaturePlots.pdf"), width = 10, height = 8)
  
  # Generate plots for each chunk
  for (chunk in feature_chunks) {
    p <- FeaturePlot(seu20, features = chunk, ncol = 2, raster = T)
    print(p)
  }
  
  # Close the PDF
  dev.off()
}

# Load Kyle's signatures and add them as a module score
kyle_signatures <- read.csv(file = "Signatures/KR_CellTypeSignatures.csv", header = T)

# Add module scores
for (n in names(kyle_signatures)) {
  print(n)
  seu20 <- AddModuleScore(object = seu20, features = kyle_signatures[n], name = n)
  colnames(seu20@meta.data) <- gsub(str_c(n, "1$"), str_c(n, "_KR_Score"), colnames(seu20@meta.data))
  }

# Define feature groups
feature_groups <- list(
  Group1 = c("Early_Erythroid_KR_Score", "Mid_Erythroid_KR_Score", "Late_Erythroid_KR_Score",
             "Monocyte_KR_Score", "cDC_KR_Score", "pDC_KR_Score", "NP_KR_Score",  "Cycling_NP_KR_Score"),
  Group2 = c("NK_KR_Score", "NKT_KR_Score", "CD56_dim_NK_KR_Score", "CD56_Bright_NK_KR_Score",
             "B_cells_KR_Score", "Pre_B_cell_KR_Score", "Plasma_Cell_KR_Score", "ProB_KR_Score"),
  Group3 = c("HSPC_KR_Score", "MPP_KR_Score", "EBM_KR_Score", "EP_KR_Score", "MkP_KR_Score","GMP_KR_Score",
             "IMP_KR_Score"), 
  Group4 = c("CD4_Naïve_KR_Score", "CD4_CM_KR_Score", "Tregs_KR_Score", "CD8_Term_Eff_KR_Score",
             "CD8_GZMK_Exh_KR_Score", "CD8_EM_KR_Score", "CD8_Naïve_KR_Score", "MAIT_KR_Score"))

# Loop through each group and save plots
for (group_name in names(feature_groups)) {
  group_features <- feature_groups[[group_name]]
  
  # Split features into smaller chunks
  feature_chunks <- split_features(group_features, features_per_page)
  
  # Open a PDF file for the group
  pdf(file = paste0("FeaturePlots/", group_name, "_Kyle_FeaturePlots.pdf"), width = 10, height = 8)
  
  # Generate plots for each chunk
  for (chunk in feature_chunks) {
    p <- FeaturePlot(seu20, features = chunk, ncol = 2, raster = T)
    print(p)
  }
  
  # Close the PDF
  dev.off()
}

# Look at Peter's published signatures (Griffin et al., 2023 - https://github.com/petervangalen/Single-cell_BPDCN/tree/main/02_Annotate)
markergenes <- read.table(file = "Signatures/Griffin_markerGenes.txt", header = T)

# Add module scores
for (n in names(markergenes)) {
  print(n)
  seu20 <- AddModuleScore(object = seu20, features = markergenes[n], name = n)
  colnames(seu20@meta.data) <- gsub(str_c(n, "1$"), str_c(n, "_PvG_Score"), colnames(seu20@meta.data))
}

feature_groups <- list(Group1 = c("HSPC_PvG_Score", "GMP_PvG_Score", "ProMono_PvG_Score",
                                  "EarlyEry_PvG_Score", "LateEry_PvG_Score"),
  Group2 = c( "Mono_PvG_Score", "ncMono_PvG_Score", "cDC_PvG_Score", "pDC_PvG_Score"), 
  Group3 = c("NK_PvG_Score", "NKT_PvG_Score", "B_PvG_Score", "PreB_PvG_Score",
             "ProB_PvG_Score", "Plasma_PvG_Score"),
  Group4 = c("CD4Naive_PvG_Score", "CD4Memory_PvG_Score", "CD8Naive_PvG_Score",
             "CD8Memory_PvG_Score", "CD8TermExh_PvG_Score", "GammaDeltaLike_PvG_Score"))

# Set the number of features per page
features_per_page <- 4

# Loop through each group and save plots
for (group_name in names(feature_groups)) {
  group_features <- feature_groups[[group_name]]
  
  # Split features into smaller chunks
  feature_chunks <- split_features(group_features, features_per_page)
  
  # Open a PDF file for the group
  pdf(file = paste0("FeaturePlots/", group_name, "_Griffin_FeaturePlots.pdf"), width = 10, height = 8)
  
  # Generate plots for each chunk
  for (chunk in feature_chunks) {
    p <- FeaturePlot(seu20, features = chunk, ncol = 2, raster = T)
    print(p)
  }
  
  # Close the PDF
  dev.off()
}


# ANNOTATE ---------------------------------------------------------------------

# Based on all the information above, add cell annotations
Idents(seu20) <- "seurat_clusters"
seu20.cluster.ids <- c("T cells", "T cells", "Monocytes", "T cells", "Late Erythroids", "T cells", "Mid Erythroids", "UD1", "B cells", "Monocytes", "Non Classical Monocytes", "T cells", "Pro Monocytes", "Late Erythroids", "Early Erythroids", "UD3", "Early Erythroids", "Pre-B", "Progenitors", "cDC", "Pro B cells", "UD1", "Early Erythroids", "Plasma cells", "Cycling T-NK cells", "pDC", "Progenitors", "Progenitors", "Pro B cells", "Late Erythroids", "UD2", "UD2", "UD2")
names(seu20.cluster.ids) <- levels(seu20)
seu20 <- RenameIdents(seu20, seu20.cluster.ids)
seu20$celltype <- Idents(seu20)

# Visualize the annotated clusters
DimPlot(seu20, reduction = "umap", repel = T, group.by = "celltype", shuffle = T, label = T, raster = T) +
  theme(aspect.ratio = 1)


# PROJECT FULL UMAP AND PREDICT CELL TYPES -------------------------------------

# Load the data from 1.1_CreateSeuratObject.R
seu <- readRDS(file = "~/250409_MergedSeuratObject.rds")

# Remove the 20% annotated cells from the full object
seu80 <- subset(seu, cells = setdiff(colnames(seu), colnames(seu20)))

# Check that the cells numbers make sense
identical(ncol(seu), ncol(seu20)+ncol(seu80))

# UMAP projection of the remaining 80% of cells (similar to https://satijalab.org/seurat/articles/integration_mapping.html#unimodal-umap-projection)
seu80 <- NormalizeData(seu80)
seu.anchors <- FindTransferAnchors(reference = seu20, query = seu80, dims = 1:18, reference.reduction = "pca")
seu80 <- IntegrateEmbeddings(anchorset = seu.anchors, reference = seu20, query = seu80, new.reduction.name = "ref.pca")
seu80 <- ProjectUMAP(query = seu80, query.reduction = "ref.pca", reference = seu20, reference.reduction = "pca", reduction.model = "umap")

# Transfer celltype annotations
predictions_df <- TransferData(anchorset = seu.anchors, refdata = seu20$celltype, dims = 1:18)
colnames(predictions_df) <- gsub("predicted.id", "celltype", colnames(predictions_df))
seu80 <- AddMetaData(seu80, metadata = select(predictions_df, celltype))

# Compare UMAPs
p2 <- DimPlot(seu20, reduction = "umap", group.by = "celltype", raster = T) + theme(aspect.ratio = 1)
p3 <- DimPlot(seu80, reduction = "ref.umap", group.by = "celltype", raster = T) + theme(aspect.ratio = 1)

pdf("2.2.2_AnnotationStage1.pdf", width = 15, height = 8)
p2 + p3
dev.off()

# Save a tibble with all cell type annotations and predictions
all_celltypes_df <- rbind(seu20@meta.data[,"celltype",drop=F], seu80@meta.data[,"celltype",drop=F])
all_celltypes_tib <- as_tibble(all_celltypes_df, rownames = "cell")
write_csv(all_celltypes_tib, file = "2.2_Stage1_celltype_annotations.csv")
gzip("2.2_Stage1_celltype_annotations.csv", overwrite = T, remove = T)

# Save a tibble with all UMAP coordinates and projections
all_coordinates_df <- rbind(seu20@reductions$umap@cell.embeddings,
                            seu80@reductions$ref.umap@cell.embeddings)
all_coordinates_tib <- as_tibble(all_coordinates_df, rownames = "cell")
write_csv(all_celltypes_tib, file = "2.2_Stage1_umap_coordinates.csv")
gzip("2.2_Stage1_umap_coordinates.csv", overwrite = T, remove = T)
