# Load the needed libraries
library(tidyverse)
library(Seurat)
library(readr)
library(stringr)

# Empty environment
rm(list=ls())

# Set working directory
setwd("~/TP53_ImmuneEscape/8_Numbat/")

# Create output folder if it doesn't exist
output_dir <- "UMI_counts"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Load the saved Seurat object
seu <- readRDS("~/250128_seurat_annotated_final.rds")

# Extract full UMI count matrix **once** before looping
seu_counts <- GetAssayData(seu, layer = "counts")

# Load the saved dataframe that contains souporcell information
final_df <- read_csv("~/final_dataset.csv")

# Extract unique patient IDs
patient_ids <- unique(final_df$patient_id)

# **Debug Step 1: Print unique patient IDs**
message("Unique patient IDs:", paste(patient_ids, collapse = ", "))

# Loop through each patient
for (patient_id in patient_ids) {
  
  message(paste("\nProcessing patient:", patient_id))  
  
  # **Fix: Properly filter `final_df` for this patient**
  df_subset <- final_df %>%
    filter(patient_id == !!patient_id, origin == "recipient") %>%  # **FIXED FILTER**
    mutate(full_barcode = paste0(orig.ident, "_", barcode))  # Ensure correct barcode format
  
  # **Debug Step 2: Check if subset is correct**
  message(paste("Total rows in `df_subset` for", patient_id, ":", nrow(df_subset)))
  if (nrow(df_subset) == 0) {
    message(paste("No matching cells found for", patient_id, "- Skipping"))
    next
  }
  
  # **Fix: Ensure only relevant barcodes are selected**
  matching_barcodes <- intersect(df_subset$full_barcode, colnames(seu_counts))
  
  # **Debug Step 3: Check number of matching barcodes**
  message(paste("Matching barcodes found for", patient_id, ":", length(matching_barcodes)))
  if (length(matching_barcodes) == 0) {
    message(paste("No matching barcodes for", patient_id, "- Skipping"))
    next
  }
  
  # **Fix: Ensure we are extracting patient-specific UMI counts**
  umi_counts <- seu_counts[, matching_barcodes, drop = FALSE]
  
  # **Debug Step 4: Check dimensions of extracted UMI counts**
  message(paste("UMI count matrix size for", patient_id, ":", dim(umi_counts)[1], "genes x", dim(umi_counts)[2], "cells"))
  
  # **Fix: If matrix is empty, skip saving**
  if (dim(umi_counts)[2] == 0) {
    message(paste("No valid UMI counts for", patient_id, "- Skipping save"))
    next
  }
  
  # Modify barcodes: Remove everything before the last underscore but keep '-1'
  colnames(umi_counts) <- sub(".*_", "", colnames(umi_counts))
  
  # Save with fast compression using `qs`
  save_path <- file.path(output_dir, paste0(patient_id, "_host_numbat_umi_counts.rds"))
  saveRDS(umi_counts, save_path)  
  
  message(paste("Saved UMI counts for", patient_id, "at", save_path))
}

# The umi counts rds file were deleted during a subsequent commit, and on 250503, also deleted from the history to reduce the size of the repository