suppressMessages(suppressWarnings({
  library(tidyverse)
  library(dplyr)
  library(Seurat)
  library(BoneMarrowMap)
  library(symphony)
  library(cowplot)
  library(MuDataSeurat)
}))

setwd("~/Documents/projects/nurefsan_integration/")

base_dir = "/Users/kr72/Documents/projects/aging_blood/prepare_for_publication/"

seu = readRDS("250416_MergedSeuratObject.rds")
samples = seu$orig.ident %>% unique()

### Map each of the samples onto BMM reference; this follows BMM tutorial ####

# Load Symphony reference (pre-downloaded)
ref = readRDS(paste0(base_dir, "input_data/BoneMarrowMap/BoneMarrow_RefMap_SymphonyRef.rds"))
# Set uwot path for UMAP projection
ref$save_uwot_path = paste0(base_dir, "input_data/BoneMarrowMap/BoneMarrow_RefMap_uwot_model.uwot")

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
bmm_annotations = bmm_annotations[,c(14:21)]

write.table(bmm_annotations, "bmm_annotations.tsv", append = F, quote = F, sep = "\t", row.names = T, col.names = T)

seu = AddMetaData(seu, bmm_annotations)

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
