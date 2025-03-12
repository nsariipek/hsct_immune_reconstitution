# Nurefsan Sariipek, modified at 250310
# Extracting recipient barcodes for running Numbat

# Load the libraries
library(tidyverse)
library(janitor)
library(readr)
library(dplyr)

# Empty environment
rm(list=ls())

setwd("~/TP53_ImmuneEscape/9_Numbat/")

# Load the saved dataframe that contains souporcell information + barcodes
final_df <- read_csv("~/final_dataset.csv")

# Get unique patient identifiers
patients <- unique(final_df$patient_id)  

# Define output directory
output_dir <- "Numbat_Barcodes/"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Loop through each patient and process their samples
for (patient in patients) {
  # Filter data for the current patient
  patient_df <- final_df %>% filter(patient_id == patient)
  
  # Extract samples by their status
  pre_transplant <- patient_df %>% filter(sample_status == "pre_transplant" & origin == "recipient")
  remission <- patient_df %>% filter(sample_status == "remission" & origin == "recipient")
  relapse <- patient_df %>% filter(sample_status == "relapse" & origin == "recipient")
  
  # Extract barcodes and remove identifier after '_'
  barcode_pre <-pre_transplant$barcode
  barcode_rem <- remission$barcode
  barcode_rel <- relapse$barcode
  
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

  
  
