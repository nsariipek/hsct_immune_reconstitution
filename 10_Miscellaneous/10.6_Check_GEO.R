# Peter van Galen, 251016
# Check that all the files for GEO submission correspond to the Seurat object names

# Load libraries
library(tidyverse)
library(readxl)
library(Seurat)

# Set working directory
repo_root <- system("git rev-parse --show-toplevel", intern = TRUE)
setwd(paste0(repo_root, "/10_Miscellaneous"))

# Clear environment variables
rm(list = ls())

# Load data
seu <- readRDS("../AuxiliaryFiles/250528_Seurat_complete.rds")
metadata <- as_tibble(seu@meta.data)

# Load files
checksums <- read_excel("10.6_GEO_checksums.xlsx", sheet = 1)
checksums

# Check fastqs ----------------------------------------------------------------

checksums_fastq_split <- checksums %>%
  filter(Folder == "fastq") %>%
  separate_wider_delim(
    File,
    delim = "_",
    names = c("patient_id", "sample_id", "library_type", "assay", "read")
  ) %>%
  select(patient_id, sample_id, library_type, assay, read)

# 30 libaries have 2 fastqs, 100 libaries have 4 fastqs
checksums_fastq_split %>%
  group_by(patient_id, sample_id, library_type, assay) %>%
  count %>%
  .$n %>%
  table

# In the fastq file list, there are 65 libraries with two assays (gex and vdj)
checksums_fastq_split %>%
  select(patient_id, sample_id, library_type, assay) %>%
  unique %>%
  group_by(patient_id, sample_id, library_type)
checksums_fastq_split %>%
  select(patient_id, sample_id, library_type, assay) %>%
  unique %>%
  .$assay %>%
  table

# The Seurat files has 65 libraries too
metasummary <- metadata %>%
  select(patient_id, sample_id, library_type) %>%
  unique %>%
  mutate(sample_id = gsub(".*_", "", sample_id))

# Finally, check if the patient, sample, and library_type data is the same
seu_ids <- select(metasummary, patient_id, sample_id, library_type) %>%
  mutate(
    patient_id = as.character(patient_id),
    library_type = as.character(library_type)
  ) %>%
  arrange(patient_id, sample_id, library_type)
fastq_ids <- select(
  checksums_fastq_split,
  patient_id,
  sample_id,
  library_type
) %>%
  unique() %>%
  arrange(patient_id, sample_id, library_type)

identical(seu_ids, fastq_ids)

# Check processed -------------------------------------------------------------
checksums_processesed_split <- checksums %>%
  filter(Folder == "processed") %>%
  separate_wider_delim(
    File,
    delim = "_",
    names = c("patient_id", "sample_id", "library_type", "file"),
    too_many = "merge"
  ) %>%
  select(patient_id, sample_id, library_type, file)

# 65 libraries have 3 file types
checksums_processesed_split %>%
  group_by(patient_id, sample_id, library_type) %>%
  count %>%
  .$n %>%
  table

# Finally, check if the patient, sample, and library_type data is the same
processed_ids <- select(
  checksums_processesed_split,
  patient_id,
  sample_id,
  library_type
) %>%
  unique() %>%
  arrange(patient_id, sample_id, library_type)

identical(seu_ids, processed_ids)


# Check processed TCR files ---------------------------------------------------

# Favorite function
cutf <- function(x, f = 1, d = "/") {
  sapply(strsplit(x, d), function(i) paste(i[f], collapse = d))
}

# This Excel sheet was generated from GEO_scHSCT_2025_PROCESSED_scTCRseq.sh
processed_scTCR <- read_excel("10.6_GEO_checksums.xlsx", sheet = 2)

# Match patient numbers from original Cell Ranger output to GEO submission file
processed_scTCR_summary <- processed_scTCR %>%
  mutate(original = cutf(original, d = "/", f = 6)) %>%
  mutate(
    final_for_geo = paste0(
      cutf(final_for_geo, d = "_|/", f = 9),
      "_",
      cutf(final_for_geo, d = "_|/", f = 11)
    )
  ) %>%
  rename(orig.ident = "original", patient_id_library_type = "final_for_geo") %>%
  arrange(orig.ident, patient_id_library_type)

# Extract orig.ident, patient ID, and library type from Seurat object
metasummary2 <- as_tibble(seu@meta.data) %>%
  select(orig.ident, patient_id, library_type) %>%
  unique %>%
  mutate(
    orig.ident = as.character(orig.ident),
    patient_id_library_type = paste0(patient_id, "_", library_type)
  ) %>%
  select(orig.ident, patient_id_library_type) %>%
  arrange(orig.ident, patient_id_library_type)

# Compare
identical(processed_scTCR_summary, metasummary2)
