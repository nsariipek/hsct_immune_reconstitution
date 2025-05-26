# Nurefsan Sariipek, 241219, updated at 250331
# Integrating Numbat Results

# Load Libraries
library(tidyverse)
library(Seurat)

# Empty environment
rm(list=ls())

# Set working directory
setwd("~/TP53_ImmuneEscape/8_Numbat/")
# For Peter:
#setwd("~/DropboxMGB/Projects/ImmuneEscapeTP53/TP53_ImmuneEscape/8_Numbat/")

# Favorite function
cutf <- function(x, f=1, d="/") sapply(strsplit(x, d), function(i) paste(i[f], collapse=d))

# Load the saved Seurat object
seu <- readRDS("../AuxiliaryFiles/250426_Seurat_annotated.rds")

# Load the saved dataframe that contains souporcell information
souporcell_assignments <- read_csv("../6_Souporcell/6.2_Souporcell_assignments.csv.gz")

# Define paths to files with Numbat data
patient_info <- list(
  P20 = "Numbat_Calls/P20_clone_post_2.tsv",
  P22 = "Numbat_Calls/P22_clone_post_2.tsv",
  P23 = "Numbat_Calls/P23_clone_post_2.tsv",
  P30 = "Numbat_Calls/P30_clone_post_2.tsv",
  P31 = "Numbat_Calls/P31_clone_post_2.tsv",
  P33 = "Numbat_Calls/P33_clone_post_2.tsv"
)

# Extract metadata from Seurat object
metadata_tib <- as_tibble(seu@meta.data, rownames = "cell")

# Loop through patients to add Numbat to the metadata
metadata_list <- list()

for (p_id in names(patient_info)) {
  #p_id <- names(patient_info)[1]
  cat("\nProcessing", p_id, "...\n")
  
  # Subset Seurat object for recipient cells of current patient
  soc_subset <- souporcell_assignments %>% filter(patient_id == p_id, origin == "recipient")
  metadata_subset_tib <- filter(metadata_tib, cell %in% soc_subset$cell)
  
  # Add a new barcode to facilitate merging with Numbat below & make sure there are no duplicates
  metadata_subset_tib$barcode <- paste0(metadata_subset_tib$patient_id, "_", cutf(metadata_subset_tib$cell, d = "_", f = 3))
  print(paste(sum(duplicated(metadata_subset_tib$barcode)), "duplicate barcodes"))

  # Load Numbat data
  clone_path <- patient_info[[p_id]]
  pt_clone <- read.table(clone_path, header = TRUE)
  # Add barcode column & select relevant columns
  pt_clone <- pt_clone %>% mutate(barcode = paste0(p_id, "_", pt_clone$cell)) %>%
    select(barcode, clone_opt, compartment_opt)
  
  # Join Seurat metadata with Numbat data & restore rownames
  meta_join <- metadata_subset_tib %>% left_join(pt_clone, by = "barcode")
  
  # Check for NAs (should be 0)
  print(paste(sum(is.na(select(meta_join, clone_opt, compartment_opt))), "NAs detected"))

  # Store
  metadata_list[[p_id]] <- meta_join
}

# Merge metadata
metadata_final <- do.call(rbind, metadata_list) %>%
  select(cell, clone_opt, compartment_opt)

# Save
write_csv(metadata_final, file = "8.4_Numbat_calls.csv")



# Check consistency with previous version
seu_combined_old <- readRDS("~/250505_numbat_combined_seurat.rds")
colnames(seu_combined_old) <- paste0(seu_combined_old$orig.ident, "_", cutf(colnames(seu_combined_old), d = "_", f = 2))
all(metadata_final$cell == colnames(seu_combined_old))
all(metadata_final$clone_opt == seu_combined_old$clone_opt)
all(metadata_final$compartment_opt == seu_combined_old$compartment_opt)

# In future scripts, add Numbat data:
numbat_calls <- read_csv("../8_Numbat/8.4_Numbat_calls.csv")
numbat_calls <- column_to_rownames(numbat_calls, var = "cell")
seu <- AddMetaData(seu, numbat_calls)