# Checking all of the samples
# Nurefsan Sariipek
# Date: "Dec 20th, 2023"
# Updated at May 27th, 2025

# Load the libraries
library(tidyverse)
library(Seurat)
library(scRepertoire)
library(ggpubr)
library(janitor)
library(rstatix)
library(Hmisc)
library(glue)
library(googleCloudStorageR) # For Terra
library(RColorBrewer)

# Empty environment
rm(list=ls())

# Set working directory 
setwd("~/TP53_ImmuneEscape/7_TCR_Diversity/")

# Setup for Peter: local. Then run line 89-117 and continue at 178
#setwd("~/DropboxMGB/Projects/ImmuneEscapeTP53/TP53_ImmuneEscape/7_TCR_Diversity")
#seu <- readRDS("../AuxiliaryFiles/250426_Seurat_annotated.rds")
#combined.sc <- readRDS("AuxiliaryFiles/combined.RDS")

# Parameters to interact with Google bucket, this part only needed for Terra
gcs_global_bucket("fc-3783b423-62ac-4c69-8c2f-98cb0ee4503b")
# Check if you can list the objects. In Terra, you may need to authenticate using gcs_auth(). In VM, this did not work - hence the alternative function on line 65.
gcs_list_objects()

# Define samples
Samples <-  c("P2446_MNC", "P25802_MNC","P25802_CD3", "P2645_MNC","P1972_MNC", "P1972_CD3", "P2220_MNC", "P2621_MNC","P2621_CD3", "P9185_MNC","P9185_CD3", "P2599_MNC","P2599_CD3", "P9596_MNC","P9596_CD3", "P25809_MNC", "P2737_MNC","P2737_CD3", "P2379_MNC","P2379_CD3", "P2434_MNC", "P2518_MNC","P2518_CD3" ,"P4618_MNC", "P6174_MNC","P6174_CD3", "P9931_MNC","P9931_CD3", "P1953_MNC","P1953_CD3", "P1677_MNC","P1677_CD3", "P1732_MNC", "P1811_MNC","P1811_CD3", "P1195_MNC", "P1285_MNC", "P1347_MNC","P1347_CD3", "P5641_MNC", "P6244_MNC","P6244_CD3", "P9355_MNC", "P1013_MNC", "P1665_MIX", "P1745_MNC", "P1817_MIX", "P2408_MNC", "P2988_MNC", "P1762_MIX", "P2698_MIX", "P2791_MNC", "P2977_MIX", "P2986_MNC", "P1671_MIX", "P2517_MIX", "P2820_MIX", "P2961_MNC", "P3000_MIX", "P1764_MIX", "P1804_MNC", "P1964_MNC", "P2332_MNC", "P2448_MNC", "P2745_MNC")

  
  # Temporary directory to save downloaded files
  tmp_dir <- "/home/rstudio/tmp"
dir.create(tmp_dir)

# Track successfully loaded samples
successful_samples <- c()

# Process each sample and assign to dynamically named variables
for (Sample in Samples) {
  print(Sample) # Log the sample being processed
  
  # Define file paths
  sample_path <- paste0(Sample, "/vdj-t/")
  file_name <- paste0(Sample, "_filtered_contig_annotations.csv")
  save_path <- file.path(tmp_dir, file_name)
  
  # Download the file and create variables
  tryCatch({
    gcs_get_object(
      object_name = paste0(sample_path, "filtered_contig_annotations.csv"),
      saveToDisk = save_path
    )
    
    # Read the downloaded file into a variable named after the sample
    assign(Sample, read.csv(save_path), envir = globalenv())
    successful_samples <- c(successful_samples, Sample) # Track success
  }, error = function(e) {
    message(glue("Error processing sample: {Sample}"))
  })
}

# Dynamically construct the contig_list for successfully loaded samples
contig_list <- mget(successful_samples, envir = globalenv())

# Combine all the samples into a single object
combined <- combineTCR(
  contig_list,
  samples = successful_samples
)

# Add variables. For scRepertoire below v2.0, replace variable.name with name
combined <- addVariable(combined, variable.name = "patient_id",
                        variables = c("P01", "P01", "P01", "P01", "P02", "P02", "P02", "P02", "P02", "P03", "P03", "P04", "P04", "P20", "P20", "P20", "P20", "P20", "P21", "P21", "P21", "P22", "P22", "P23", "P23", "P23", "P23", "P23", "P23", "P23", "P30", "P30", "P30", "P30", "P30", "P31", "P31", "P31", "P31", "P32", "P32", "P32", "P33", "P33", "P07", "P05", "P08", "P09", "P06", "P10", "P11", "P12", "P13", "P14", "P15", "P16", "P17", "P18", "P19", "P24", "P25", "P28", "P27", "P29", "P26"))


combined <- addVariable(combined, variable.name = "timepoint",
                        variables = c("pre-tx", "remission", "remission", "remission", "pre-tx", "pre-tx", "remission", "remission", "remission", "remission", "remission", "remission", "remission", "pre-tx", "pre-tx", "remission", "relapse", "relapse", "pre-tx", "pre-tx", "remission", "remission", "remission", "pre-tx", "remission", "remission", "remission", "remission", "relapse", "relapse", "pre-tx", "pre-tx", "remission", "relapse", "relapse", "pre-tx", "remission", "relapse", "relapse", "pre-tx", "remission", "remission", "pre-tx", "remission", "remission", "remission", "remission", "remission", "remission", "remission", "remission", "remission", "remission", "remission", "remission", "remission", "remission", "remission", "remission", "remission", "remission", "remission", "remission", "remission", "remission"))

# Optional: merge data if the same sample was analyzed as both MNC and sorted T cells
combined2 <- do.call(rbind, combined)
combined <- split(combined2, f = combined2$patient_id)

# load the T cells
Tcells <- readRDS("~/250428_Tcells.rds")

# Turn to a dataframe and keep only needed variables
meta = Tcells@meta.data 
meta$barcode <- rownames(meta)
meta %>%
  select(barcode, celltype, cohort, sample_status, orig.ident, sample_id, patient_id,timepoint, library_type,TP53_status)
meta = meta %>% drop_na()

  