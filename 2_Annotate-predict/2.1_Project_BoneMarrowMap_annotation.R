# Nurefsan is trying to run the BM annotation using Ksenia's script, 250424
# Load the libraries
  library(tidyverse)
  library(dplyr)
  library(Seurat)
  library(BoneMarrowMap)
  library(symphony)
  library(RColorBrewer)
  library(patchwork)
  
# For running the BM library first time run this :  
  # ## install dependencies that are not on CRAN
  # if(!require(BiocManager, quietly = TRUE)) install.packages("BiocManager")
  # BiocManager::install(c("AUCell", "doMC", "BiocNeighbors"))
  # if(!require(devtools, quietly = TRUE)) install.packages("devtools")
  # devtools::install_github("jaredhuling/jcolors")
  # ## install BoneMarrowMap package
  # devtools::install_github('andygxzeng/BoneMarrowMap', force = TRUE)

# Set working directory
setwd("~/TP53_ImmuneEscape/2_Annotate-predict/")

# Load the object that Peter advised to use
seu = readRDS("~/250417_MergedSeuratObject.rds")

# Set directory to store projection reference files
projection_path = '~/'

# Download Bone Marrow Reference - 344 Mb
curl::curl_download('https://bonemarrowmap.s3.us-east-2.amazonaws.com/BoneMarrowMap_SymphonyReference.rds', 
                    destfile = paste0(projection_path, 'BoneMarrowMap_SymphonyReference.rds'))
# Download uwot model file - 221 Mb
curl::curl_download('https://bonemarrowmap.s3.us-east-2.amazonaws.com/BoneMarrow_RefMap_uwot_model.uwot', 
                    destfile = paste0(projection_path, 'BoneMarrowMap_uwot_model.uwot'))

# Load Symphony reference
ref <- readRDS(paste0(projection_path, 'BoneMarrowMap_SymphonyReference.rds'))
# Set uwot path for UMAP projection
ref$save_uwot_path <- paste0(projection_path, 'BoneMarrowMap_uwot_model.uwot')

# ReferenceSeuratObj <- create_ReferenceObject(ref)
# DimPlot(ReferenceSeuratObj, reduction = 'umap', group.by = 'CellType_Annotation_formatted', 
#         raster=FALSE, label=TRUE, label.size = 4) + NoAxes()

### Map each of the samples onto BMM reference; this follows BMM tutorial ####
samples <- unique(seu$orig.ident)

for (i in 1:length(samples)) {
  
  query = subset(seu, orig.ident==samples[i])
  
  cat(sprintf("Processing sample %s, %d/%d, %d cells", samples[i], i, length(samples), ncol(query)), "\n")
  
  # Map query dataset using Symphony 
  query = map_Query(
    exp_query = query@assays$RNA$counts,
    metadata_query = query@meta.data,
    ref_obj = ref
  )
  
  # Run QC based on mapping error score, flag cells with mapping error >= 2.5 MADs above median
  query = query %>% calculate_MappingError(., reference = ref, MAD_threshold = 2.5) 
  
  # Predict Hematopoietic Cell Types by KNN classification
  query = predict_CellTypes(
    query_obj = query, 
    ref_obj = ref, 
    initial_label = 'initial_CellType', # celltype assignments before filtering on mapping QC
    final_label = 'predicted_CellType'  # celltype assignments with map QC failing cells assigned as NA
  )
  query = AddMetaData(query, query@reductions$umap@cell.embeddings)
  
  write.table(query@meta.data, paste0(samples[i], ".bmm_mapped.tsv"), append = F, quote = F, sep = "\t", row.names = T, col.names = T)
  
  gc()
  
}

annotation_files = list.files(".", pattern = ".*bmm_mapped.tsv", recursive = TRUE, full.names = T)

bmm_annotations = data.frame()

for (i in 1:length(annotation_files)) {
  t = read.table(annotation_files[[i]], header = T, sep = "\t", row.names = 1)
  bmm_annotations = rbind(bmm_annotations, t)
}

bmm_annotations = bmm_annotations[colnames(seu), ]
all(colnames(seu) == rownames(bmm_annotations))

# Select the needed columns and saved them as a tsv
bmm_annotations = bmm_annotations[,c(16:23)]

write.table(bmm_annotations, "bmm_annotations.tsv", append = F, quote = F, sep = "\t", row.names = T, col.names = T)

# Add into your seurat object
seu = AddMetaData(seu, bmm_annotations)

# Take out the UMAP coordinates
umap_coords = as.matrix(seu@meta.data[, c("umap_1", "umap_2")])

# Ensure row names match the Seurat object cell names
rownames(umap_coords) = colnames(seu)

seu[["umap_bmm"]] = CreateDimReducObject(
  embeddings = umap_coords,
  key = "umap_bmm_", 
  assay = DefaultAssay(seu)
)

DimPlot(seu, group.by = "predicted_CellType", label = T) & scale_color_manual(values = c(brewer.pal("Paired", n=12), brewer.pal(n=11, "Spectral"), brewer.pal(n=9, "Set1"), "black", "purple", "grey", "cyan", "gold", "green", "violet", brewer.pal("Paired", n=12), brewer.pal("Paired", n=12))) & NoLegend()

DimPlot(seu, group.by = "predicted_CellType", label = F, split.by = "patient_id", ncol=6) & scale_color_manual(values = c(brewer.pal("Paired", n=12), brewer.pal(n=11, "Spectral"), brewer.pal(n=9, "Set1"), "black", "purple", "grey", "cyan", "gold", "green", "violet", brewer.pal("Paired", n=12), brewer.pal("Paired", n=12))) & NoLegend()

saveRDS(seu, "~/250424_seu_bm_annotated.rds")

