# Nurefsan Sariipek, 241219, updated at 250331
# Integrating Numbat Results

# Load Libraries
library(readr)
library(Seurat)
library(tidyverse)
library(RColorBrewer)
library(ggsci)

# Empty environment
rm(list=ls())

# Set working directory
setwd("~/TP53_ImmuneEscape/8_Numbat/")
# For Peter:
#setwd("~/DropboxMGB/Projects/ImmuneEscapeTP53/TP53_ImmuneEscape/8_Numbat/")

# Favorite function
cutf <- function(x, f=1, d="/") sapply(strsplit(x, d), function(i) paste(i[f], collapse=d))

# Load the saved Seurat objects
seu <- readRDS("../AuxiliaryFiles/250426_Seurat_annotated.rds")

# Extract UMAP coordinates from the Seurat object
umap_coords <- Embeddings(seu, reduction = "umap_bmm")

# Rename the columns for clarity
colnames(umap_coords) <- c("UMAP_1", "UMAP_2")

# Add them to metadata
seu@meta.data <- cbind(seu@meta.data, umap_coords)

# Load the saved dataframe that contains souporcell information
souporcell_assignments <- read_csv("../6_Souporcell/6.2_Souporcell_assignments.csv.gz")

# Helper function
make_unique <- function(x) {
  make.unique(x, sep = "__")
}

# Define Numbat clone paths
patient_info <- list(
  P20 = "Numbat_Calls/P20_clone_post_2.tsv",
  P22 = "Numbat_Calls/P22_clone_post_2.tsv",
  P23 = "Numbat_Calls/P23_clone_post_2.tsv",
  P30 = "Numbat_Calls/P30_clone_post_2.tsv",
  P31 = "Numbat_Calls/P31_clone_post_2.tsv",
  P33 = "Numbat_Calls/P33_clone_post_2.tsv"
)

# Step 1 — Patch souporcell_assignments barcodes
# --------------------------
souporcell_assignments$barcode <- paste0(souporcell_assignments$patient_id, "_",
                                         cutf(souporcell_assignments$cell, d = "_", f = 3))
souporcell_assignments$barcode <- make_unique(souporcell_assignments$barcode)

# --------------------------
# Step 2 — Patch Seurat colnames
# --------------------------
barcode_df <- as.data.frame(seu@meta.data)
barcode_df$old_barcode <- rownames(barcode_df)
barcode_df$patient_id <- barcode_df$patient_id %||% barcode_df$Sample %||% NA

if (any(is.na(barcode_df$patient_id))) {
  stop("❌ Cannot find 'patient_id' or 'Sample' in Seurat metadata.")
}

new_barcodes <- paste0(barcode_df$patient_id, "_", sub(".*_", "", barcode_df$old_barcode))
new_barcodes <- make_unique(new_barcodes)

colnames(seu) <- new_barcodes
rownames(seu@meta.data) <- new_barcodes

# --------------------------
# Step 3 — Loop through patients
# --------------------------
seu_list <- list()

for (patient_id in names(patient_info)) {
  cat("\n▶ Processing", patient_id, "...\n")
  
  clone_path <- patient_info[[patient_id]]
  if (!file.exists(clone_path)) {
    warning(paste("⚠ Numbat file not found for", patient_id))
    next
  }
  
  # Subset souporcell
  soc_subset <- souporcell_assignments %>%
    filter(patient_id == !!patient_id, origin == "recipient")
  
  matching_cells <- soc_subset$barcode
  available_cells <- intersect(matching_cells, colnames(seu))
  
  if (length(available_cells) == 0) {
    warning(paste("⚠ No matching barcodes in Seurat for", patient_id))
    next
  }
  
  seu_subset <- subset(seu, cells = available_cells)
  
  # Load Numbat
  pt_clone <- read.table(clone_path, header = TRUE)
  pt_clone$barcode <- paste0(patient_id, "_", pt_clone$cell)
  pt_clone$barcode <- make_unique(pt_clone$barcode)
  
  pt_select <- pt_clone %>%
    select(barcode, clone_opt, compartment_opt)
  
  # Join metadata
  meta_to_add <- tibble(barcode = colnames(seu_subset)) %>%
    left_join(pt_select, by = "barcode")
  
  if (any(is.na(meta_to_add$clone_opt))) {
    warning(paste("⚠ Some barcodes did not match between Seurat and Numbat for", patient_id))
    next
  }
  
  seu_subset <- AddMetaData(seu_subset, metadata = meta_to_add %>% column_to_rownames("barcode"))
  seu_subset$patient_id <- patient_id
  
  # Store
  seu_list[[patient_id]] <- seu_subset
}

# --------------------------
# Step 4 — Merge all
# --------------------------
if (length(seu_list) == 0) {
  stop("❌ No patients were successfully processed.")
} else if (length(seu_list) == 1) {
  seu_combined <- seu_list[[1]]
  cat("\n✅ Only one patient processed — no merging needed.\n")
} else {
  cat("\n✅ Merging all patients...\n")
  seu_combined <- merge(x = seu_list[[1]], y = seu_list[-1], project = "CombinedPatients")
}

# --------------------------
# Step 5 — Save
# --------------------------
saveRDS(seu_combined, "~/250505_numbat_combined_seurat.rds")
cat("\n🎉 All done! Saved as 'combined_patients_seurat.rds'\n")


