# Nurefsan Sariipek and Peter van Galen, 250518
# To run Souporcell on merged bam files, we are merging the barcodes from each patients here, and we are removing the barcodes that are present in multiple samples
# See https://github.com/wheaton5/souporcell for souporcell explanation

# Load the libraries
library(tidyverse)
library(Seurat)

# Start with a clean slate
rm(list = ls())

# Set the working directory (for Nurefsan)
setwd("~/hsct_immune_reconstitution/5_Souporcell/")
# For Peter:
setwd(
  "~/DropboxMGB/Projects/ImmuneEscapeTP53/hsct_immune_reconstitution/08_Souporcell/"
)

# Load Seurat object
seu <- readRDS("../AuxiliaryFiles/250426_Seurat_annotated.rds")

# Extract metadata
df1 <- seu@meta.data
df1 <- df1 %>% rownames_to_column("barcode")

# Select patient IDs
patient_ids <- unique(df1$patient_id)

# Function to process each patient
process_patient <- function(patient) {
  cat("\nProcessing:", patient, "...\n") # Print progress

  df_filtered <- df1 %>%
    filter(patient_id == patient) %>%
    mutate(barcode = str_extract(barcode, "[^_]+$")) %>%
    count(barcode, name = "Freq")

  # Count duplicates & unique barcodes
  num_duplicates <- sum(df_filtered$Freq > 1)
  num_unique <- sum(df_filtered$Freq == 1)
  total_excluded <- num_duplicates

  # Print progress
  cat(" - Duplicate barcodes:", num_duplicates, "\n")
  cat(" - Unique barcodes:", num_unique, "\n")
  cat(" - Total excluded:", total_excluded, "\n")

  # Save only unique barcodes
  df_filtered %>%
    filter(Freq == 1) %>%
    select(barcode) %>%
    write_tsv(
      file = paste0("barcodes/", patient, ".barcodes.tsv"),
      col_names = FALSE
    )

  cat("✅ Done processing", patient, "\n")
}

# Apply function to all patients & track progress
walk(patient_ids, process_patient)
