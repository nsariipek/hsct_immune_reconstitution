# Nurefsan Sariipek, created at 230405, updated at 250113
# Load the libraries
library(tidyverse)
library(Seurat)
#devtools::install_github('immunogenomics/presto')

# Start with a clean slate
rm(list=ls())

# Load the saved Seurat object from the last step that contains only 20% of the cells
seu <- readRDS("~/250113_SplittedSeuratObject.rds")

# Run FindAllMarkers to see the most expressed genes
seu_markers <- FindAllMarkers(seu, min.pct = .3, logfc.threshold = .3)

# And save as CVS since it will be easier to work it on
seu_markers_tib <- as_tibble(seu_markers)
write.csv(seu_markers_tib, file = "~/250113_marker_genes.csv")

#First visualize the current clusters and save as PDF, this will guide you
umap_plot <- DimPlot(seu, reduction = "umap", group.by = "seurat_clusters", shuffle = T, label = T) + theme(aspect.ratio = 1,legend.position = "none")
pdf("UMAP_Cluster_Labels.pdf", width = 10, height = 8)  # Adjust width and height as needed
print(umap_plot)
dev.off()

#General Features
FeaturePlot(seu, features = c("CD34","MPO", "CD14", "MS4A1", "CD3G","CD8B"))

# Some features to start with 
#T cell Features
FeaturePlot(seu, features = c("CD8A", "TCF7", "TOX", "HAVCR2", "CXCR3",
                              "SLAMF6", "CD3E", "CD4", "SELL", "CD44", "PDCD1",
                              "FOXP3", "GZMB", "GZMK", "LAG3", "CD101",
                              "CXCR5", "KLRG1", "IFNG", "TNF")) 

#HSC markers
FeaturePlot(seu, features = c("VIM", "FLT3", "CD34", "ITGAL", "THY1",
                              "PTPRC", "KIT", "SLAMF1", "MME", "SLAMF2", "MPO")) 
#B cell markers
FeaturePlot(seu, features = c("CD19", "MS4A1", "PDCD1LG2", "NT5E", "FCER2", "SDC1",
                              "PAX5", "TCF3", "CD80", "SPIB", "BCL6")) 

# Since there is still no universal way to annotate cells, we have used previously annotated single cell analysis from our lab to help us with the annotations
# Unpublish data from discovery cohort

discoverygenes <- read.table(file = "~/TP53_ImmuneEscape/2_Annotate/Markers/Cohort1_MarkerGenes.txt", 
                             sep = "\t", 
                             header = TRUE,  # Matches col.names = TRUE in write.table
                             quote = "",     # Matches quote = FALSE in write.table
                             stringsAsFactors = FALSE)  # Optional, prevents automatic factor conversion

for (n in names(discoverygenes)) {
  print(n)
  #n <- "HSPC"
  seu <- AddModuleScore(object = seu, features = discoverygenes[n], name = n)
  colnames(seu@meta.data) <- gsub(str_c(n, "1$"), str_c(n, "_Score"), colnames(seu@meta.data))
}

# Visualize
FeaturePlot(seu, features = c("B.cells_Score", "CD4.Memory_Score", "CD4.Naïve_Score", "CD56.Bright.NK.cells_Score", "CD56.Dim.NK.cells_Score", "CD8.Central.Memory_Score", "CD8.Effector.Memory_Score", "CD8.Naïve_Score", "CD8.Terminally.Exhausted_Score", "cDC_Score", "Doublets_Score", "Early.Erythroids_Score", "HSPCs_Score", "Late.Erythroids_Score", "Mid.Erythroids_Score", "Monocytes_Score", "NK.T.cells_Score", "Non.Classical.Monocytes_Score", "pDC_Score", "Plasma.Cells_Score", "Pre.B.cells_Score", "Pro.B.cells_Score", "Pro.Monocytes_Score", "Treg_Score", "Undetermined_Score", "γδ.T.lymphocytes_Score"
))


# Function to split features into chunks
split_features <- function(features, chunk_size) {
  split(features, ceiling(seq_along(features) / chunk_size))
}

# Define your specific feature groups
feature_groups <- list(
  Group1 = c("CD4.Memory_Score", "CD4.Naïve_Score", "CD56.Bright.NK.cells_Score", 
             "CD56.Dim.NK.cells_Score", "Treg_Score", "γδ.T.lymphocytes_Score", 
             "CD8.Central.Memory_Score", "CD8.Effector.Memory_Score", "CD8.Naïve_Score", 
             "NK.T.cells_Score", "CD8.Terminally.Exhausted_Score"),
  Group2 = c("Doublets_Score", "Early.Erythroids_Score", "pDC_Score", 
             "Late.Erythroids_Score", "Mid.Erythroids_Score", "cDC_Score"),
  Group3 = c("Monocytes_Score", "Non.Classical.Monocytes_Score", 
             "Pro.Monocytes_Score", "Undetermined_Score", "HSPCs_Score"),
  Group4 = c("Pre.B.cells_Score", "Pro.B.cells_Score", 
             "B.cells_Score", "Plasma.Cells_Score")
)

# Set the number of features per page
features_per_page <- 4

# Loop through each group and save plots
for (group_name in names(feature_groups)) {
  group_features <- feature_groups[[group_name]]
  
  # Split features into smaller chunks
  feature_chunks <- split_features(group_features, features_per_page)
  
  # Open a PDF file for the group
  pdf(file = paste0(group_name, "_FeaturePlots.pdf"), width = 10, height = 8)
  
  # Generate plots for each chunk
  for (chunk in feature_chunks) {
    p <- FeaturePlot(seu, features = chunk, ncol = 2)
    print(p)
  }
  
  # Close the PDF
  dev.off()
}


# Upload Kyle's signatures and add them as a module score
Kylecells <- read.csv(file = "~/KR_CellTypeSignatures.csv", header = T)

# Add Module Scores 
for (n in names(Kylecells)) {
  print(n)
  #n <- "HSPC"
seu_diet <- AddModuleScore(object = seu_diet, features = Kylecells[n], name = n)
  colnames(seu_diet@meta.data) <- gsub(str_c(n, "1$"), str_c(n, "_Score"), colnames(seu_diet@meta.data))}

# Visualize
FeaturePlot(seu_diet, features = c("NK_Score","Monocyte_Score","B_cells_Score","Early_Erythroid_Score", "Mid_Erythroid_Score", "Late_Erythroid_Score", "cDC_Score", "pDC_Score", "Plasma_Cell_Score", "Pre_B_cell_Score", "Cycling_NP_Score", "NP_Score", "ProB_Score", "GMP_Score", "HSPC_Score", "EP_Score", "IMP_Score", "MkP_Score", "MPP_Score", "EBM_Score", "CD4_Naïve_Score", "CD56_dim_NK_Score", "CD8_Term_Eff_Score", "CD8_GZMK_Exh_Score", "CD8_EM_Score", "CD8_Naïve_Score", "NKT_Score", "CD4_CM_Score", "MAIT_Score", "Tregs_Score", "CD56_Bright_NK_Score"))

# Upload Erica's signatures and add them as a module score
# Published dataset see ....
markergenes <- read.table(file = "~/markerGenes.txt", header = T)

for (n in names(markergenes)) {
  print(n)
  #n <- "HSPC"
  seu_diet <- AddModuleScore(object = seu_diet, features = markergenes[n], name = n)
  colnames(seu_diet@meta.data) <- gsub(str_c(n, "1$"), str_c(n, "_Score"), colnames(seu_diet@meta.data))
}

# Visualize
FeaturePlot(seu_diet, features = c("HSPC_Score",	"EarlyEry_Score",	"LateEry_Score",	
                              "GMP_Score",	"ProMono_Score",	"Mono_Score",	"ncMono_Score",	
                              "cDC_Score",	"pDC_Score",	"ProB_Score",	"PreB_Score",	
                              "B_Score",	"Plasma_Score",	"CD4Naive_Score",	"CD4Memory_Score",
                              "CD8Naive_Score",	"CD8Memory_Score",	"CD8TermExh_Score",
                              "GammaDeltaLike_Score",	"NKT_Score",	"NK_Score"))




###### Rename clusters #####
View(seu_diet@meta.data)

Idents(seu_diet) = "seurat_clusters"

# Rename Clusters
seu.cluster.ids <- c("T cells","T cells","T cells","Monocytes","Blasts","T cells","T cells","Late Erythroids","Monocytes","Mid Erythroids", "B cells","T cells","Pro Monocytes", "Mid Erythroids","Mid Erythroids","Late Erythroids","Blasts","Early Erythroids","Non Classical Monocytes", "B cells", "Monocytes", "Monocytes", "cDC","Early Erythroids","Pre B cells", "HSPCs", "Unidentified", "Plasma Cells", "Blasts", "Blasts", "Pro B cells", "Late Erythroids", "pDC", "Late Erythroids", "Doublets")

names(seu.cluster.ids) <- levels(seu_diet)
seu_diet <- RenameIdents(seu_diet, seu.cluster.ids)
seu_diet@meta.data$celltype = Idents(seu_diet)

# See the levels
levels(seu_diet$celltype)
levels(seu_diet)
View(seu_diet@meta.data)

# Visualize the annotated clusters
mycolors <- distinctColorPalette(k = 34)
pie(rep(1, 34), col = mycolors) 
DimPlot(seu_diet_merged, reduction = "umap", repel = T, group.by = "celltype", label = T) + theme(aspect.ratio = 1)
