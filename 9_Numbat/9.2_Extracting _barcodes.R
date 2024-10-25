# Nurefsan Sariipek, 241018, 
# Extracting receipient barcodes for running Numbat
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

# Load the saved dataframe that contaions souporcell information + barcodes
combined_df <- read_csv(paste0(my_wd,file = "/results/cohort1-2_souporcell.csv"))

#Select the sample you want the get the barcodes and select only the receipient cells

pt5 <- combined_df %>%
      subset(subset = orig.ident %in% c("2737_MNC","25809_MNC","9596_MNC"))  %>% 
      subset(subset = assignment =="host")
        

pt5_rel <- combined_df %>%
           subset(subset = orig.ident =="2737_MNC") %>% 
           subset(subset = assignment =="host")


pt5_rem <- combined_df %>%
          subset(subset = orig.ident =="25809_MNC") %>% 
          subset(subset = assignment =="host")

pt5_pre <-combined_df %>%
  subset(subset = orig.ident =="9596_MNC") %>% 
  subset(subset = assignment =="host")



#Select the barcodes
barcode_pt5_pre <- pt5_pre$cell

#Remove the sample identifier after _
barcode_pt5_pre = gsub("_.*", "", barcode_pt5_pre)
head(barcode_pt5_pre)

#export as tsv file since souporcell requires that
write_tsv((as.data.frame(barcode_pt5_pre)), file = "~/pt5_pre_host_barcodes.tsv", col_names = FALSE)


