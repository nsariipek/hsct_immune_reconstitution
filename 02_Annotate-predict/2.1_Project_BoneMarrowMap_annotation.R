# Nurefsan Sarripek, 250424
# Run BM annotation using Ksenia's script
# The final version was run by Peter on a Google VM on 250426

# Load the libraries
library(tidyverse)
library(Seurat)
library(BoneMarrowMap)
library(symphony)
library(RColorBrewer)
library(patchwork)
library(ggsci)

# For running the BM library first time:
# ## install dependencies that are not on CRAN
# if(!require(BiocManager, quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install(c("AUCell", "doMC", "BiocNeighbors"))
# if(!require(devtools, quietly = TRUE)) install.packages("devtools")
# devtools::install_github("jaredhuling/jcolors")
# ## install BoneMarrowMap package
# devtools::install_github('andygxzeng/BoneMarrowMap')

# Set working directory
setwd("~/hsct_immune_reconstitution/2_Annotate-predict/")

# Clear environment
rm(list = ls())

# Load the object generated in 1_Seurat/1.1_CreateSeuratObject.R
seu <- readRDS("../AuxiliaryFiles/250417_MergedSeuratObject.rds")


### RUN BONE MARROW MAP CELL TYPE PREDICTION ###

# Set directory to store projection reference files
projection_path = '~/'

# Download Bone Marrow Reference - 344 Mb
curl::curl_download(
  'https://bonemarrowmap.s3.us-east-2.amazonaws.com/BoneMarrowMap_SymphonyReference.rds',
  destfile = paste0(projection_path, 'BoneMarrowMap_SymphonyReference.rds')
)
# Download uwot model file - 221 Mb
curl::curl_download(
  'https://bonemarrowmap.s3.us-east-2.amazonaws.com/BoneMarrow_RefMap_uwot_model.uwot',
  destfile = paste0(projection_path, 'BoneMarrowMap_uwot_model.uwot')
)

# Load Symphony reference
ref <- readRDS(paste0(projection_path, 'BoneMarrowMap_SymphonyReference.rds'))
# Set uwot path for UMAP projection
ref$save_uwot_path <- paste0(projection_path, 'BoneMarrowMap_uwot_model.uwot')

# Visualize reference (optional)
ReferenceSeuratObj <- create_ReferenceObject(ref)
DimPlot(
  ReferenceSeuratObj,
  reduction = 'umap',
  group.by = 'CellType_Annotation_formatted',
  raster = FALSE,
  label = TRUE,
  label.size = 4
) +
  NoAxes() +
  theme(aspect.ratio = 1, legend.position = "none")

### Map each of the samples onto BMM reference; this follows BMM tutorial ####
samples <- unique(seu$orig.ident)

for (i in 1:length(samples)) {
  query = subset(seu, orig.ident == samples[i])

  cat(
    sprintf(
      "Processing sample %s, %d/%d, %d cells",
      samples[i],
      i,
      length(samples),
      ncol(query)
    ),
    "\n"
  )

  # Map query dataset using Symphony
  query = map_Query(
    exp_query = query@assays$RNA$counts,
    metadata_query = query@meta.data,
    ref_obj = ref
  )

  # Run QC based on mapping error score, flag cells with mapping error >= 2.5 MADs above median
  query = query %>%
    calculate_MappingError(., reference = ref, MAD_threshold = 2.5)

  # Predict hematopoietic cell types by KNN classification
  query = predict_CellTypes(
    query_obj = query,
    ref_obj = ref,
    initial_label = 'initial_CellType', # celltype assignments before filtering on mapping QC
    final_label = 'predicted_CellType' # celltype assignments with map QC failing cells assigned as NA
  )
  query = AddMetaData(query, query@reductions$umap@cell.embeddings)

  # Add pseudotime value prediction by KNN
  query <- predict_Pseudotime(
    query_obj = query,
    ref_obj = ref,
    initial_label = 'initial_Pseudotime', # pseudotime assignments before filtering on mapping QC
    final_label = 'predicted_Pseudotime' # pseudotime assignments with map QC failing cells assigned as NA
  )

  write.table(
    query@meta.data,
    paste0(samples[i], ".bmm_mapped.tsv"),
    append = F,
    quote = F,
    sep = "\t",
    row.names = T,
    col.names = T
  )

  gc()
}


### ADD ANNOTATIONS TO SEURAT OBJECT ###

# Merge the annotation files and order as Seurat object
annotation_files = list.files(
  ".",
  pattern = ".*bmm_mapped.tsv",
  recursive = TRUE,
  full.names = T
)

bmm_annotations = data.frame()

for (i in 1:length(annotation_files)) {
  t = read.table(annotation_files[[i]], header = T, sep = "\t", row.names = 1)
  bmm_annotations = rbind(bmm_annotations, t)
}

bmm_annotations = bmm_annotations[colnames(seu), ]
all(colnames(seu) == rownames(bmm_annotations))

# Select the needed columns and saved them as a tsv
bmm_annotations = bmm_annotations[, c(14:25)]

write.table(
  bmm_annotations,
  "../AuxiliaryFiles/250426_BMM_annotations.tsv",
  append = F,
  quote = F,
  sep = "\t",
  row.names = T,
  col.names = T
)

# Add predicted cell type annotations into the Seurat object
seu = AddMetaData(
  seu,
  bmm_annotations[, c(
    "mapping_error_QC",
    "predicted_CellType",
    "predicted_CellType_Broad",
    "predicted_Pseudotime"
  )]
)

# Take out the UMAP coordinates
umap_coords = as.matrix(bmm_annotations[, c("umap_1", "umap_2")])

# Ensure row names match the Seurat object cell names
all(rownames(umap_coords) == colnames(seu))

seu[["umap_bmm"]] = CreateDimReducObject(
  embeddings = umap_coords,
  key = "umapBMM_",
  assay = DefaultAssay(seu)
)

# Check how many cells did not pass QC
as_tibble(seu@meta.data) %>%
  group_by(patient_id, mapping_error_QC) %>%
  count() %>%
  group_by(patient_id) %>%
  mutate(proportion = n / sum(n)) %>%
  ggplot(aes(x = patient_id, y = proportion, fill = mapping_error_QC)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#D53E4F", "#4DAF4A")) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    panel.grid = element_blank(),
    aspect.ratio = 0.5
  )
ggsave("2.1.1_QC_pass-fail.pdf", width = 8, height = 4)


### DEFINE CELL TYPE ORDER AND COLOR ###

# Create a color vector
#mycol <-  c(brewer.pal("Paired", n=12), brewer.pal(n=11, "Spectral"), brewer.pal(n=9, "Set1"), brewer.pal(n=8, "Set2"), brewer.pal(n=12, "Set3"), brewer.pal(n=8, "Accent"), pal_igv("default")(51))
#mycol <- mycol[!duplicated(mycol)]

# Create cell type vector
#celltype_colors_tib <- tibble(celltype = c(unique(c(seu$predicted_CellType, #seu$predicted_CellType_Broad)), rep(NA, 36)), color = mycol)
#celltype_colors_tib <- celltype_colors_tib %>%
#  mutate(granularity = case_when(
#    celltype %in% seu$predicted_CellType & celltype %in% seu$predicted_CellType_Broad ~ "Both",
#    celltype %in% seu$predicted_CellType ~ "Granular",
#    celltype %in% seu$predicted_CellType_Broad ~ "Broad",
#    .default = NA
#  ))
# The following file was edited and moved to the root directory
#write_tsv(celltype_colors_tib, "celltype_colors.txt")

# Load colors
celltype_colors_df <- read.table(
  "../celltype_colors.txt",
  sep = "\t",
  header = T,
  stringsAsFactors = F,
  comment.char = ""
)
celltype_colors <- setNames(
  celltype_colors_df$color,
  celltype_colors_df$celltype
)


### FINALIZE ANNOTATION ###

# We will move forward with the broad cell types with the following exceptions
add_granular_celltypes <- c(
  "CD14 Mono",
  "CD16 Mono",
  "Immature B",
  "Mature B",
  "CD4 Naive",
  "CD4 Central Memory",
  "CD4 Effector Memory",
  "CD4 Regulatory",
  "CD8 Naive",
  "CD8 Central Memory",
  "CD8 Effector Memory 1",
  "CD8 Effector Memory 2",
  "CD8 Tissue Resident Memory",
  "T Proliferating",
  "NK CD56high",
  "NK Proliferating"
)

# Create final celltype column
seu@meta.data <- seu@meta.data %>%
  mutate(
    celltype = if_else(
      predicted_CellType %in% add_granular_celltypes,
      predicted_CellType, # use granular if in list
      predicted_CellType_Broad # otherwise use broad
    )
  )

# Then order them logically
seu$predicted_CellType <- factor(
  seu$predicted_CellType,
  levels = intersect(celltype_colors_df$celltype, seu$predicted_CellType)
)
seu$predicted_CellType_Broad <- factor(
  seu$predicted_CellType_Broad,
  levels = intersect(celltype_colors_df$celltype, seu$predicted_CellType_Broad)
)
seu$celltype <- factor(
  seu$celltype,
  levels = intersect(celltype_colors_df$celltype, seu$celltype)
)

# Check
#celltype_counts_tib <- as_tibble(seu@meta.data) %>% group_by(predicted_CellType,
#  predicted_CellType_Broad, celltype) %>% count() %>% view

### PLOT UMAPS ###

# Local for Peter. This only works for the "celltype" group below
setwd(
  "~/DropboxMGB/Projects/ImmuneEscapeTP53/hsct_immune_reconstitution/2_Annotate-predict/"
)
#seu <- readRDS("../AuxiliaryFiles/250426_Seurat_annotated.rds")
#annotated_cells <- colnames(seu)[!is.na(seu$celltype)]
#seu_pass <- subset(seu, cells = annotated_cells)

# Subset for cells with a good annotation
seu_pass <- subset(seu, mapping_error_QC == "Pass")

# Plot projected cells
DimPlot(
  seu_pass,
  group.by = "predicted_CellType",
  label = T,
  raster = T,
  raster.dpi = c(1536, 1536),
  pt.size = 3
) +
  scale_color_manual(values = celltype_colors) +
  theme(aspect.ratio = 1)
ggsave("2.1.2_Granular_annotation.pdf", width = 16, height = 8)

DimPlot(
  seu_pass,
  group.by = "celltype",
  label = T,
  raster = T,
  raster.dpi = c(1536, 1536),
  pt.size = 3
) +
  scale_color_manual(values = celltype_colors) +
  theme(aspect.ratio = 1)
ggsave("2.1.3_Final_annotation.pdf", width = 12, height = 8)

DimPlot(
  seu_pass,
  group.by = "celltype",
  label = F,
  split.by = "patient_id",
  ncol = 6
) +
  scale_color_manual(values = celltype_colors) +
  theme(
    panel.border = element_rect(color = "black"),
    aspect.ratio = 1,
    title = element_blank()
  )
ggsave("2.1.4_Per_patient_annotation.pdf", width = 20, height = 16)


### SAVE AND CLEAN UP ###

# Remove unnecessary information (it's still in 250426_BMM_annotations.tsv if needed)
seu$mapping_error_QC <- NULL
seu$predicted_CellType <- NULL
seu$predicted_CellType_Broad <- NULL

# Save finalized Seurat object
saveRDS(seu, "../AuxiliaryFiles/250426_Seurat_annotated.rds")
#system("gcloud storage cp ../AuxiliaryFiles/250426_Seurat_annotated.rds gs://fc-3783b423-62ac-4c69-8c2f-98cb0ee4503b")

# Clean up
#unlink("~/BoneMarrowMap_SymphonyReference.rds")
#unlink("~/BoneMarrowMap_uwot_model.uwot")
#unlink(annotation_files)
