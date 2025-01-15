# Nurefsan Sariipek, created at 230405, updated at 250113
# Load the libraries
library(tidyverse)
library(Seurat)

# Start with a clean slate
rm(list=ls())

# Load the saved Seurat object from the last step that containes only 20% of the cells
seu <- readRDS("~/250113_SplittedSeuratObject.rds")

# Run FindAllMarkers to see the most expressed genes
seu_markers <- FindAllMarkers(seu, min.pct = .3, logfc.threshold = .3)

# Save as tibble, as Peter requested
write.tsv(as_tibble(seu_markers), file = "~/250113_marker_genes.tsv")

# And save as CVS since it will be easier to work it on
seu_markers_tib <- as_tibble(seu_markers)
write.csv(seu_markers_tib, file = "~/250113_marker_genes.csv")

# Load marker genes that you have saved before
seu_markers <- read.csv(file = "~/250113_marker_genes.csv")

# To visualize and help with annotating clusters run the lines below
#General Features
FeaturePlot(seu_diet, features = c("CD34","MPO", "CD14", "MS4A1", "CD3G","CD8B"))

# Feature Plots
#T cell Features
FeaturePlot(seu_diet, features = c("CD8A", "TCF7", "TOX", "HAVCR2", "CXCR3",
                              "SLAMF6", "CD3E", "CD4", "SELL", "CD44", "PDCD1",
                              "FOXP3", "GZMB", "GZMK", "LAG3", "CD101",
                              "CXCR5", "KLRG1", "IFNG", "TNF")) 
#Myeloid and B cell Features
FeaturePlot(seu_diet, features = c("CCR2", "CD14", "CD33", "CD34", "THY1", "IL3RA",
                              "CX3CR1", "MME", "PTPRC", "ITGAX", "CD80", "CD19",
                              "ITGAL", "XCR1", "LY6C", "CSF1R", "ADGRE1")) 
#HSC markers
FeaturePlot(seu_diet, features = c("VIM", "FLT3", "CD34", "ITGAL", "THY1",
                              "PTPRC", "KIT", "SLAMF1", "MME", "SLAMF2", "MPO")) 
#B cell markers
FeaturePlot(seu_diet, features = c("CD19", "MS4A1", "PDCD1LG2", "NT5E", "FCER2", "SDC1",
                              "PAX5", "TCF3", "CD80", "SPIB", "BCL6")) 


# Since there is still no universal way to annotate cells we have used previously annotated single cell anaylysis from our lab to help us with the annotations
# Unpublish data 
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


#Also we haved used several of published datasets, for annotation purposes on Terra which @Nurefsan will try to add those on here in the future.

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
