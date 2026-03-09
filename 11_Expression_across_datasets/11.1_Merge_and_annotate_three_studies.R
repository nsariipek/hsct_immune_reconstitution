# Peter van Galen, 260110
# Merge and annotate three datasets (the goal is to visualize gene expression in the next script)

library(tidyverse)
library(Seurat)
library(BoneMarrowMap)
library(symphony)
library(ggsci)

# Set working directory
repo_root <- system("git rev-parse --show-toplevel", intern = TRUE)
setwd(paste0(repo_root, "/AnalysisPeter"))

# Clear environment variables & load favorite function
rm(list = ls())
cutf <- function(x, f = 1, d = "/") {
  sapply(strsplit(x, d), function(i) paste(i[f], collapse = d))
}

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

# Load data from three different projects
# Van Galen et al., Cell 2019. Object is available from Figshare (https://doi.org/10.6084/m9.figshare.30581066)
seu_aml_og <- readRDS("../AuxiliaryFiles/Seurat_AML.rds")
# Griffin et al., Nature 2023. Object is available from Figshare (https://doi.org/10.6084/m9.figshare.31467496)
seu_bpdcn_og <- readRDS("../AuxiliaryFiles/BM_Seurat_Final.rds")
# Sariipek et al. (this work). Object is available from Figshare (https://doi.org/10.6084/m9.figshare.30213781)
seu_hsct_og <- readRDS("../AuxiliaryFiles/250528_Seurat_complete.rds")

# Common genes
common_genes <- intersect(
  rownames(seu_aml_og),
  intersect(rownames(seu_bpdcn_og), rownames(seu_hsct_og))
)

# Get data from 2019 single-cell AML paper ------------------------------------

seu_aml <- subset(seu_aml_og, features = common_genes)
seu_aml <- CreateSeuratObject(
  counts = LayerData(seu_aml, assay = "RNA", layer = "counts"),
  meta.data = seu_aml@meta.data
)

# Filter for healthy control (BM) or malignant AML cells at Dx.
bm_ids <- seu_aml@meta.data |> filter(grepl("^BM", orig.ident)) |> rownames()
aml_ids <- seu_aml@meta.data |>
  filter(grepl("^AML.*D0", orig.ident), PredictionRefined == "malignant") |>
  rownames()
seu_aml <- subset(seu_aml, cells = c(bm_ids, aml_ids))

# Only keep necessary metadata columns
seu_aml@meta.data <- data.frame(
  orig.ident = seu_aml@meta.data$orig.ident,
  project = "aml_2019",
  tech = "Seq-well",
  donor_id = cutf(as.character(seu_aml@meta.data$orig.ident), d = "\\."),
  replicate = seu_aml@meta.data$orig.ident,
  cell_origin = case_when(
    grepl("BM", seu_aml@meta.data$orig.ident) ~ "Healthy",
    grepl("AML", seu_aml@meta.data$orig.ident) ~ "AML"
  ),
  select(seu_aml@meta.data, PredictionRefined),
  CellType1 = seu_aml@meta.data$CellType
) |>
  # Add leading 0 to BM IDs so they differ from the 2023 study
  mutate(donor_id = gsub("BM", "BM0", donor_id))

# Get healthy bone marrow donors from BPDCN project ---------------------------

seu_bpdcn <- subset(seu_bpdcn_og, features = common_genes)
seu_bpdcn <- CreateSeuratObject(
  counts = LayerData(seu_bpdcn, assay = "RNA", layer = "counts"),
  meta.data = seu_bpdcn@meta.data
)

# Only keep necessary metadata columns
seu_bpdcn@meta.data <- data.frame(
  orig.ident = seu_bpdcn@meta.data$orig.ident,
  project = "healthy_donors_2023",
  # fmt: skip
  tech = gsub("TenX", "10X_3prime",
    gsub("SW", "Seq-well", seu_bpdcn@meta.data$tech)
  ),
  donor_id = cutf(seu_bpdcn@meta.data$replicate, d = "\\."),
  replicate = seu_bpdcn@meta.data$replicate,
  cell_origin = "Healthy",
  CellType2 = seu_bpdcn@meta.data$CellType,
  row.names = colnames(seu_bpdcn)
) |>
  # Add five to donor ID so they differ from the 2019 study
  mutate(donor_id = sprintf("BM%02d", as.numeric(gsub("BM", "", donor_id)) + 5))

# Get data from 2026 HSCT paper -----------------------------------------------

seu_hsct <- subset(seu_hsct_og, features = common_genes)
seu_hsct <- CreateSeuratObject(
  counts = LayerData(seu_hsct, assay = "RNA", layer = "counts"),
  meta.data = seu_hsct@meta.data
)

# Filter for cells with celltype and souporcell assignments
exclude_ids <- seu_hsct_og@meta.data |>
  filter(is.na(celltype) | is.na(souporcell_origin)) |>
  rownames()
seu_hsct <- subset(seu_hsct, cells = setdiff(colnames(seu_hsct), exclude_ids))

# Also exclude early relapse cohort for better agreement with the HSCT paper
seu_hsct <- subset(seu_hsct, cohort_detail != "1-Early-relapse")

# Only keep necessary metadata columns
seu_hsct@meta.data <- data.frame(
  orig.ident = seu_hsct@meta.data$orig.ident,
  project = "hsct_2026",
  tech = "10X_5prime",
  donor_id = seu_hsct@meta.data$patient_id,
  select(
    seu_hsct@meta.data,
    sample_id,
    library_type,
    cohort,
    cohort_detail,
    timepoint,
    sample_status,
  ),
  cell_origin = case_when(
    seu_hsct@meta.data$souporcell_origin == "recipient" ~ "Recipient",
    seu_hsct@meta.data$souporcell_origin == "donor" ~ "Donor"
  ),
  numbat_compartment = seu_hsct@meta.data$numbat_compartment,
  CellType3 = seu_hsct@meta.data$celltype
)


# Merge -----------------------------------------------------------------------

# Merge three datasets together
seu_merge <- merge(seu_aml, merge(seu_bpdcn, seu_hsct))
seu_merge <- JoinLayers(seu_merge)
seu_merge <- NormalizeData(seu_merge)


# Annotate --------------------------------------------------------------------

# This and following sections are similar to the procedures in 02_Annotate-predict/2.1_Project_BoneMarrowMap_annotation.R
# This section can be skipped if it was run before and the annotation files are present

# Download BoneMarrowMap_SymphonyReference.rds (344 Mb) and BoneMarrowMap_uwot_model.uwot (221 Mg)
curl::curl_download(
  'https://bonemarrowmap.s3.us-east-2.amazonaws.com/BoneMarrowMap_SymphonyReference.rds',
  destfile = '../AuxiliaryFiles/BoneMarrowMap_SymphonyReference.rds'
)
curl::curl_download(
  'https://bonemarrowmap.s3.us-east-2.amazonaws.com/BoneMarrow_RefMap_uwot_model.uwot',
  destfile = '../AuxiliaryFiles/BoneMarrowMap_uwot_model.uwot'
)

# Load Symphony reference
ref <- readRDS('../AuxiliaryFiles/BoneMarrowMap_SymphonyReference.rds')
# Set uwot path for UMAP projection
ref$save_uwot_path <- '../AuxiliaryFiles/BoneMarrowMap_uwot_model.uwot'

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

# Map each of the samples onto BMM reference; this follows BMM tutorial
samples <- unique(seu_merge$donor_id)

for (i in 1:length(samples)) {
  query = subset(seu_merge, donor_id == samples[i])

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


# Add annotations to Seurat object --------------------------------------------

# Merge the annotation files and order as Seurat object
annotation_files <- list.files(
  ".",
  pattern = ".*bmm_mapped.tsv",
  recursive = T,
  full.names = T
)

bmm_annotations <- data.frame()

for (i in 1:length(annotation_files)) {
  t = read.table(annotation_files[[i]], header = T, sep = "\t", row.names = 1)
  bmm_annotations = rbind(bmm_annotations, t)
}

bmm_annotations <- bmm_annotations[colnames(seu_merge), ]
all(colnames(seu_merge) == rownames(bmm_annotations))

# Select the needed columns and saved them as a tsv
bmm_annotations <- bmm_annotations[, c(20:31)]

write_tsv(
  rownames_to_column(bmm_annotations, "cell"),
  "../AuxiliaryFiles/260110_BMM_annotations_ThreeStudies.tsv"
)

# Add predicted cell type annotations into the Seurat object
seu_merge <- AddMetaData(
  seu_merge,
  bmm_annotations[, c(
    "mapping_error_QC",
    "predicted_CellType",
    "predicted_CellType_Broad",
    "predicted_Pseudotime"
  )]
)

# Take out the UMAP coordinates
umap_coords = as.matrix(bmm_annotations[, c(
  "umapprojected_1",
  "umapprojected_2"
)])

# Ensure row names still match the Seurat object cell names
all(rownames(umap_coords) == colnames(seu_merge))

seu_merge[["umap_bmm"]] = CreateDimReducObject(
  embeddings = umap_coords,
  key = "umapBMM_",
  assay = DefaultAssay(seu_merge)
)

# Check how many cells did not pass QC
as_tibble(seu_merge@meta.data) %>%
  group_by(donor_id, mapping_error_QC) %>%
  count() %>%
  group_by(donor_id) %>%
  mutate(proportion = n / sum(n)) %>%
  ggplot(aes(x = donor_id, y = proportion, fill = mapping_error_QC)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#D53E4F", "#4DAF4A")) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    panel.grid = element_blank(),
    aspect.ratio = 0.5
  )
ggsave("11.1.1_QC_pass-fail.pdf", width = 8, height = 4)


# Finalize annotation ---------------------------------------------------------

# As in 02_Annotate-predict/2.1_Project_BoneMarrowMap_annotation.R, we will move forward with broad cell types with the following exceptions
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
seu_merge@meta.data <- seu_merge@meta.data %>%
  mutate(
    celltype = if_else(
      predicted_CellType %in% add_granular_celltypes,
      predicted_CellType, # use granular if in list
      predicted_CellType_Broad # otherwise use broad
    )
  )

# Then order them logically
seu_merge$predicted_CellType <- factor(
  seu_merge$predicted_CellType,
  levels = intersect(celltype_colors_df$celltype, seu_merge$predicted_CellType)
)
seu_merge$predicted_CellType_Broad <- factor(
  seu_merge$predicted_CellType_Broad,
  levels = intersect(
    celltype_colors_df$celltype,
    seu_merge$predicted_CellType_Broad
  )
)
seu_merge$celltype <- factor(
  seu_merge$celltype,
  levels = intersect(celltype_colors_df$celltype, seu_merge$celltype)
)

# Check
#celltype_counts_tib <- as_tibble(seu_merge@meta.data) %>% group_by(predicted_CellType,
#  predicted_CellType_Broad, celltype) %>% count() %>% view

# Plot UMAPs ------------------------------------------------------------------

# Subset for cells with a good annotation
seu_pass <- subset(seu_merge, mapping_error_QC == "Pass")

# Plot projected cells
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
ggsave("11.1.2_Final_annotation.pdf", width = 12, height = 8)

DimPlot(
  seu_pass,
  group.by = "celltype",
  pt.size = 3,
  label = F,
  split.by = "donor_id",
  ncol = 6
) +
  scale_color_manual(values = celltype_colors) +
  theme(
    panel.border = element_rect(color = "black"),
    aspect.ratio = 1,
    title = element_blank()
  )
ggsave("11.1.3_Per_patient_annotation.pdf", width = 28, height = 28)


# Save and clean up -----------------------------------------------------------

# Remove unnecessary information (it's still in 260110_BMM_annotations_ThreeStudies.tsv if needed)
seu_pass$mapping_error_QC <- NULL
seu_pass$predicted_CellType <- NULL
seu_pass$predicted_CellType_Broad <- NULL

# Save finalized Seurat object
saveRDS(seu_pass, "../AuxiliaryFiles/260110_Seurat_ThreeStudies.rds")

# Clean up
unlink("../AuxiliaryFiles/BoneMarrowMap_SymphonyReference.rds")
unlink("../AuxiliaryFiles/BoneMarrowMap_uwot_model.uwot")
unlink(annotation_files)

# Run in terminal:
#gcloud storage cp AuxiliaryFiles/260110_Seurat_ThreeStudies.rds gs://fc-3783b423-62ac-4c69-8c2f-98cb0ee4503b
