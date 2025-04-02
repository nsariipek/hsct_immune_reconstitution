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
setwd("~/TP53_ImmuneEscape/9_Numbat/")

# Load the saved Seurat objects
seu <- readRDS("~/250128_seurat_annotated_final.rds")


# Load the saved dataframe that contains souporcell information
final_Df <- read_csv("~/final_dataset.csv")


# Helper function

make_unique <- function(x) {
  make.unique(x, sep = "__")
}

# Define Numbat clone paths

patient_info <- list(
  P01 = "Numbat_Calls/P01_clone_post_2.tsv",
  P05 = "Numbat_Calls/P05_clone_post_2.tsv",
  P07 = "Numbat_Calls/P07_clone_post_2.tsv",
  P08 = "Numbat_Calls/P08_clone_post_2.tsv",
  P09 = "Numbat_Calls/P09_clone_post_2.tsv",
  P10 = "Numbat_Calls/P10_clone_post_2.tsv",
  P12 = "Numbat_Calls/P12_clone_post_2.tsv",
  P13 = "Numbat_Calls/P13_clone_post_2.tsv",
  P14 = "Numbat_Calls/P14_clone_post_2.tsv",
  P17 = "Numbat_Calls/P17_clone_post_2.tsv",
  P18 = "Numbat_Calls/P18_clone_post_2.tsv",
  P23 = "Numbat_Calls/P23_clone_post_2.tsv",
  P24 = "Numbat_Calls/P24_clone_post_2.tsv",
  P25 = "Numbat_Calls/P25_clone_post_2.tsv",
  P30 = "Numbat_Calls/P30_clone_post_2.tsv"
)

# Step 1 — Patch final_Df barcodes
# --------------------------
final_Df$barcode <- paste0(final_Df$patient_id, "_", final_Df$barcode)
final_Df$barcode <- make_unique(final_Df$barcode)

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
  soc_subset <- final_Df %>%
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
saveRDS(seu_combined, "~/numbat_combined_seurat.rds")
cat("\n🎉 All done! Saved as 'combined_patients_seurat.rds'\n")


