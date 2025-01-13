# Nurefsan Sariipek, 241219
# Combine Numbat seurat object for uMAP visulazation


# Load Libraries
library(readr)
library(numbat)
library(Seurat)
library(tidyverse)
library(RColorBrewer)
library(ggsci)

# Empty environment
rm(list=ls())

# Set working directory
# For Nurefsan
setwd("/Users/dz855/Dropbox (Partners HealthCare)/ImmuneEscapeTP53/TP53_ImmuneEscape/9_Numbat")

numbat_seu <- readRDS("/Users/dz855/Partners HealthCare Dropbox/Nurefsan Sariipek/ImmuneEscapeTP53/TP53_ImmuneEscape/9_Numbat/numbat_seu.RDS")

# Load subsetted seurat object for each patient 
seu_pt5 <- readRDS("/Users/dz855/Partners HealthCare Dropbox/Nurefsan Sariipek/ImmuneEscapeTP53/TP53_ImmuneEscape/9_Numbat/pt5.RDS")
seu_pt8 <- readRDS("/Users/dz855/Partners HealthCare Dropbox/Nurefsan Sariipek/ImmuneEscapeTP53/TP53_ImmuneEscape/9_Numbat/pt8.RDS")
seu_pt9 <- readRDS("/Users/dz855/Partners HealthCare Dropbox/Nurefsan Sariipek/ImmuneEscapeTP53/TP53_ImmuneEscape/9_Numbat/pt9.RDS")
seu_pt10 <- readRDS("/Users/dz855/Partners HealthCare Dropbox/Nurefsan Sariipek/ImmuneEscapeTP53/TP53_ImmuneEscape/9_Numbat/pt10.RDS")
seu_pt12 <- readRDS("/Users/dz855/Partners HealthCare Dropbox/Nurefsan Sariipek/ImmuneEscapeTP53/TP53_ImmuneEscape/9_Numbat/pt12.RDS")

# Merge
numbat_seu <- merge(seu_pt5, y = list(seu_pt8, seu_pt9, seu_pt10, seu_pt12),add.cell.ids = c("pt5", "pt8", "pt9","pt10","pt12"))

# Add prefixes to the UMAP embeddings for each object
umap_pt5 <- Embeddings(seu_pt5, "umap")
rownames(umap_pt5) <- paste0("pt5_", rownames(umap_pt5))

umap_pt8 <- Embeddings(seu_pt8, "umap")
rownames(umap_pt8) <- paste0("pt8_", rownames(umap_pt8))

umap_pt9 <- Embeddings(seu_pt9, "umap")
rownames(umap_pt9) <- paste0("pt9_", rownames(umap_pt9))

umap_pt10 <- Embeddings(seu_pt10, "umap")
rownames(umap_pt10) <- paste0("pt10_", rownames(umap_pt10))

umap_pt12 <- Embeddings(seu_pt12, "umap")
rownames(umap_pt12) <- paste0("pt12_", rownames(umap_pt12))

# Now, combine the UMAP embeddings from all the original objects:
  
# Combine all UMAP embeddings
combined_umap <- rbind(umap_pt5, umap_pt8, umap_pt9, umap_pt10, umap_pt12)

# Add the combined UMAP embeddings to the merged Seurat object:
numbat_seu[["umap"]] <- CreateDimReducObject(
  embeddings = combined_umap,
  key = "UMAP_",
  assay = DefaultAssay(numbat_seu)
)

# Plot UMAP
pdf("UMAP_plot.pdf", width = 7, height = 7)  # Adjust width and height as needed
DimPlot(numbat_seu, reduction = "umap", group.by = "compartment_opt",cols = c("#228B22","#B22222"), split.by = "status") +
  theme(aspect.ratio = 1)
dev.off()


# Generate the UMAP plot with celltype annotations
cluster_sizes <- table(numbat_seu@meta.data$celltype)
print(cluster_sizes)
threshold <- 500
clusters_to_label <- names(cluster_sizes[cluster_sizes > threshold])
print(clusters_to_label)

umap_plot <- DimPlot(
  numbat_seu, 
  reduction = "umap", 
  group.by = "celltype", 
  split.by = "compartment_opt") + 
  theme(aspect.ratio = 1, legend.position = "none")

# Add labels on top of clusters
umap_plot <- LabelClusters(
  plot = umap_plot, 
  id = "celltype",  
  clusters = clusters_to_label,
  repel = TRUE      )

# Display the plot
umap_plot

ggsave("UMAP_with_larger_split_titles.pdf", plot = umap_plot, width = 12, height = 8, dpi = 300)

# Save Final merged numbat seurat object
# saveRDS(numbat_seu,"/Users/dz855/Partners HealthCare Dropbox/Nurefsan Sariipek/ImmuneEscapeTP53/TP53_ImmuneEscape/9_Numbat/numbat_seu.RDS")
