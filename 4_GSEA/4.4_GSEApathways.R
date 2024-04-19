# Visualizing the GSEA pathways as a module score on my dataset

# Add a GMT file as a signature score to a Seurat object in R
# Load necessary libraries:
library(Seurat)
library(qusage)

# Clean the environment
rm(list=ls())

# Load your Seurat object:
my_wd <- "/Users/dz855/Dropbox (Partners HealthCare)/ImmuneEscapeTP53/"
# Load the seurat object
seu1 <- readRDS(paste0(my_wd, "AnalysisNurefsan/RDS files/seu_diet_merged.rds"))

#Load the T cells
#to use with assignments
#seu <- readRDS(paste0(my_wd, "AnalysisNurefsan/RDS files/assignment_seu.rds"))

# Only select the post-Tx 3-6 months samples
seu1 <- subset(x=seu1, subset = Sample %in% c("P01_1Rem", "P01_2Rem", "P02_1Rem", "P04_1Rem",  "P05_1Rem", "P06_1Rem", "P07_1Rem", "P08_1Rem"))
# Select NK cells
seu <- subset(x=seu, subset = celltype %in% c("CD56 Bright NK cells","CD56 Dim NK cells"))

# Select only T cells and NK cells
seu <- subset(x=seu1, subset = celltype %in% c("CD8 Effector", "CD8 Terminally Exhausted", "NK T cells","CD8 Memory","γδ T lymphocytes","CD4 Memory","Treg", "CD4 Naïve","CD8 Naïve","CD56 Dim NK cells", "CD56 Bright NK cells"))


seu <- subset(x=seu1, subset = celltype %in% c("CD56 Dim NK cells"))

# Select donor cells
#seu <- subset(x=seu, subset = assignment =="donor")

# Load your gene set from the GMT file
gmt_file <- "/Users/dz855/Dropbox (Partners HealthCare)/ImmuneEscapeTP53/AnalysisNurefsan/DGE/signatures/GOMF_G_PROTEIN_COUPLED_RECEPTOR_ACTIVITY.v2023.2.Hs.gmt"
  #"/Users/dz855/Dropbox (Partners HealthCare)/ImmuneEscapeTP53/AnalysisNurefsan/DGE/signatures/REACTOME_EUKARYOTIC_TRANSLATION_ELONGATION.v2023.2.Hs.gmt"
 # "/Users/dz855/Dropbox (Partners HealthCare)/ImmuneEscapeTP53/AnalysisNurefsan/DGE/signatures/GOMF_LIGAND_GATED_CHANNEL_ACTIVITY.v2023.2.Hs.gmt"
  #"/Users/dz855/Dropbox (Partners HealthCare)/ImmuneEscapeTP53/AnalysisNurefsan/DGE/signatures/IL15_UP.V1_DN.v2023.2.Hs.gmt"

  
# Convert to a list 
gene_sets <- read.gmt(gmt_file)

# Or add a genelist from a paper
cd4 <- read.csv(paste0(my_wd,"AnalysisNurefsan/DGE/signatures/cd4.csv"), sep = "\t", header = T)
cd8 <- read.csv(paste0(my_wd,"AnalysisNurefsan/DGE/signatures/cd8.csv"), sep = "\t", header = T)
neoantigen <- read.csv(paste0(my_wd,"AnalysisNurefsan/DGE/signatures/neoantigen.csv"), sep = "\t", header = T)
# Add to the seurat object
neoantigen<- AddModuleScore(object = Tcells_combined, 
                     features = neoantigen,
                     name = "neoantigen",
                     assay = "RNA",
                     search = T)

# View 
View(neoantigen@meta.data)

# Plot it
# Change the identity
neoantigen <- SetIdent(neoantigen, value = "Sample")

# Visualize it 
FeaturePlot(seu, features = "ADGRG1", split.by = "cohort")

VlnPlot(seu, features = "ADGRG1", split.by = "cohort", sort = "increasing")

# Define an order of cluster identities
#my_levels <- c(4,3,2,1)

#seu@ident <- factor(x = object@ident, levels = my_levels)

#### Troubleshooting ####

# 1. Overlapping genes from GSEA and DGE and plot these genes indivudially 
# Turn GSEA gene sets to data frame and rename the column
t1 <- as.data.frame(gene_sets)
colnames(t1) <- "gene"
# Load the DE results (paired--The same one I have used for the GSEA)
test1 <- read_csv(paste0(my_wd, "AnalysisNurefsan/DGE/only_Tcells-nonsig/NKcellsonly.pseudobulk_DE_res.csv"))
# Filter out the ones have pvalue> 0.05
# Set thresholds
# padj_cutoff = 0.05
# test2 <- filter(test1, padj < padj_cutoff)

# Find the overlapping genes and turn into df
common_values <- intersect(t1$gene, test1$gene)
t2 <- as.data.frame(common_values)
test4 <- test1 %>% 
         filter(gene %in% common_values) %>% 
         arrange(desc(log2FoldChange))

t5 <- print(head(test4$gene))

gene_sets <-as.list(t5)

# Turn to a list to add a module score to code in previous section or  
gene_sets <- as.list(t2)


# 2. Take the first(sorted by logfoldchane) 10 or genes from the GSEA pathway and visualize them instead of the all genes either as a signature or individually
# Here are my thoughts 
# How do I know which ones are the first 10 ?
# Should I look into DGE again, maybe ?
# But isn't that the same thing that i did on the number 1 ?


# Load the GSEA 
test3 <- read_tsv(paste0(my_wd, "AnalysisNurefsan/DGE/only_Tcells-nonsig/GSEA_NKcells_result.tsv"))





