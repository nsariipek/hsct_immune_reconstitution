# Nurefsan Sariipek, modified at 250310
# Extracting recipient barcodes for running Numbat

# Load the libraries
library(tidyverse)
library(janitor)
library(readr)
library(dplyr)

setwd("~/TP53_ImmuneEscape/8_Numbat/")

# Empty environment
rm(list=ls())

# Favorite function
cutf <- function(x, f=1, d="/") sapply(strsplit(x, d), function(i) paste(i[f], collapse=d))

# Load the saved dataframe that contains souporcell information + barcodes
souporcell_assignments <- read_csv("../6_Souporcell/6.2_Souporcell_assignments.csv.gz")

# Get unique patient identifiers
patients <- unique(souporcell_assignments$patient_id)

# Define output directory
dir.create("Numbat_Barcodes/")

# Loop through each patient and process their samples
for (patient in patients) {
  # Filter data for the current patient
  patient_df <- souporcell_assignments %>% filter(patient_id == patient)
  
  # Extract samples by their status
  pre_transplant <- patient_df %>% filter(sample_status == "pre-transplant" & origin == "recipient")
  remission <- patient_df %>% filter(sample_status == "remission" & origin == "recipient")
  relapse <- patient_df %>% filter(sample_status == "relapse" & origin == "recipient")
  
  # Extract barcodes 
  barcode_pre <- cutf(pre_transplant$cell, d = "_", f = 3)
  barcode_rem <- cutf(remission$cell, d = "_", f = 3)
  barcode_rel <- cutf(relapse$cell, d = "_", f = 3)
  
  # Convert to tibble
  barcode_pre <- as_tibble(barcode_pre)
  barcode_rem <- as_tibble(barcode_rem)
  barcode_rel <- as_tibble(barcode_rel)
  
  # Save barcodes as TSV files
  write_tsv(as.data.frame(barcode_pre), file = paste0(output_dir, patient, "_pre_recipient_barcodes.tsv"), col_names = FALSE)
  write_tsv(as.data.frame(barcode_rem), file = paste0(output_dir, patient, "_rem_recipient_barcodes.tsv"), col_names = FALSE)
  write_tsv(as.data.frame(barcode_rel), file = paste0(output_dir, patient, "_rel_recipient_barcodes.tsv"), col_names = FALSE)
  
  print(paste("Saved barcode files for patient:", patient))
}

# On 250502, we deleted the Numbat_Barcodes folder, since the Patient IDs are now outdated and we would need to regenerate the barcode files to rerun Numbat, which are not planning to do