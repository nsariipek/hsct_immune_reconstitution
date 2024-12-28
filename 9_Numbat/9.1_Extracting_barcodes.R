# Nurefsan Sariipek, 241018 
# Extracting recipient barcodes for running Numbat
# Adjust code for each patient separately, subsetting recipient cells each time

# Load the libraries
library(tidyverse)
library(janitor)
library(Seurat)
library(readr)
library(dplyr)
library(numbat)
library(fossil)
library(readxl)

# Empty environment
rm(list=ls())

# For Nurefsan:
my_wd <- "/Users/dz855/Dropbox (Partners HealthCare)/ImmuneEscapeTP53/TP53_ImmuneEscape/5_Souporcell/"
# For Peter
my_wd <- "~/DropboxMGB/Projects/ImmuneEscapeTP53/TP53_ImmuneEscape/5_Souporcell/"

# Load the saved dataframe that contaions souporcell information + barcodes
combined_df <- read_csv(paste0(my_wd,file = "/results/cohort1-2_souporcell.csv"))

# Load the following for the cohort 3
combined_df <- read_csv(paste0(my_wd,file = "/results/cohort3_souporcell.csv"))

#Select the sample you want the get the barcodes and select only the recipient cells
pt10_pre <- combined_df %>%
           subset(subset = orig.ident =="1195_MNC") %>% 
           subset(subset = assignment =="host")

pt10_rem <- combined_df %>%
          subset(subset = orig.ident =="1285_MNC") %>% 
          subset(subset = assignment =="host")

pt10_rel <- combined_df %>%
  subset(subset = orig.ident =="1347_MNC") %>% 
  subset(subset = assignment =="host")

# pt11_rel <- combined_df %>%
#            subset(subset = orig.ident =="1347_MNC") %>% 
#            subset(subset = assignment =="host")

# Select the barcodes
barcode_pt10_rel_host <- pt10_rel$cell

# Remove the sample identifier after _
barcode_pt10_pre_host  = gsub("_.*", "", barcode_pt10_pre_host)
head(barcode_pt10_pre_host)
barcode_pt10_pre_host <- as_tibble(barcode_pt10_pre_host)

# Export as tsv file since Numbat requires that
write_tsv((as.data.frame(barcode_pt10_rel_host)), file = "~/pt10_rel_host_barcodes.tsv", col_names = FALSE)


