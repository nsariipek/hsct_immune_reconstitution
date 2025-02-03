# Nurefsan Sariipek, 230920+ Peter van Galen, 230925
# Updated at 250203
# Since we have merged samples from different time points for each patient, we have to merge the barcodes as well to be able to run the Souporcell.
# See https://github.com/wheaton5/souporcell to understand how it works

# Load the libraries
library(tidyverse)
library(Seurat)

# Start with a clean slate
rm(list=ls())

# Set the working directory 
setwd("~/TP53_ImmuneEscape/5_Souporcell/")

# Load Seurat object
seu <- readRDS("~/250128_seurat_annotated_final.rds")

# Since we merged the bam files from each patient, we are merging the barcodes from each patient as well, and we are removing the barcodes that are present more than 2 sample.
df1 <- seu@meta.data
df1 = df1 %>% rownames_to_column("barcode")

# Select patient IDs 
patient_ids <- unique(df1$patient_id)

# Function to process each patient
process_patient <- function(patient) {
  cat("\nProcessing:", patient, "...\n")  # Print progress
  
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
    write_tsv(file = paste0("barcodes/", patient, ".barcodes.tsv"), col_names = FALSE)
  
  cat("✅ Done processing", patient, "\n")
}

# Apply function to all patients & track progress
walk(patient_ids, process_patient)


#Old code-delete at some point
# 
# # Select all the samples belonging to patient 9.
# P9 <-  subset(x = df1, subset = patient_id=="P09")
# 
# #Select the barcodes
# barcode9 <- P9$barcode
# # Remove everything before the last underscore
# barcode9 <- gsub(".*_", "", barcode9)
# 
# # Print the result
# print(barcode9)
# 
# 
# # Exclude the doublets that across the different samples
# barcode9_dups <- as.data.frame(table(barcode9)) %>% filter(Freq > 1) %>% pull(barcode9)
# barcode9_unique <- as.data.frame(table(barcode9)) %>% filter(Freq == 1) %>% pull(barcode9)
# 
# # Check the number of rows to see how many dropped 
# length(barcode9_dups) # 238 barcodes occur more than once
# length(barcode9_unique) # 21,001 barcodes occur once
# length(barcode9) - length(barcode9_unique) # 479 cells will be excluded
# 
# #export as tsv file since souporcell requires that
# write_tsv((as.data.frame(barcode9_unique)), file = "barcodes/pt9.barcodes.tsv", col_names = FALSE)
# 
# # Overlap between all combinations
# P9_barcodes_ls <- P9 %>% select(barcode, orig.ident) %>%
#   mutate(barcode = gsub("-.*", "", barcode)) %>%
#   select(orig.ident, barcode) %>%
#   group_by(orig.ident) %>% 
#   pivot_wider(names_from = orig.ident, values_from = barcode, values_fn = list)


# # Addition from Peter
# 
# # Make a list of every sample compared to every other sample
# overlap_list <- lapply(P9_barcodes_ls, function(x) lapply(P9_barcodes_ls, function(y) sum(duplicated(c(x[[1]], y[[1]])))))
# overlap_list
# 
# # Generate a dataframe
# do.call(cbind, overlap_list)
