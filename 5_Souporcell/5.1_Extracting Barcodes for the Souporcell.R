# Nurefsan Sariipek, 230920
# Peter van Galen, 230925
# Since we have merged samples from different time points for each patient, we have to merge the barcodes as well to be able to run the Souporcell.
# See https://github.com/wheaton5/souporcell to understand how it works

# Load the libraries
library(tidyverse)
library(Seurat)

# Start with a clean slate
rm(list=ls())

# Set the working directory 
# For Nurefsan:
my_wd <- "/Users/dz855/Dropbox (Partners HealthCare)/ImmuneEscapeTP53/"

# Load Seurat object
seu_diet_merged <- readRDS(paste0(my_wd,file ="RDS files/seu_diet_merged.rds"))

# Since we merged the bam files from each patient, we are merging the barcodes from each patient as well, and we are removing the barcodes that are present more than 2 sample.
df1 <- seu_diet_merged@meta.data
df1 = df1 %>% rownames_to_column("barcode")

# Select all the samples belonging to patient 9.
P9 <-  subset(x = df1, subset = orig.ident %in% c("1677_MNC", "1732_MNC","1811_MNC","1677_CD3","1811_CD3"))

#Select the barcodes
barcode9 <- P9$barcode

#Remove the sample identifier after _
barcode9 = gsub("_.*", "", barcode9)

# Exclude the doublets that across the different samples
barcode9_dups <- as.data.frame(table(barcode9)) %>% filter(Freq > 1) %>% pull(barcode9)
barcode9_unique <- as.data.frame(table(barcode9)) %>% filter(Freq == 1) %>% pull(barcode9)

# Check the number of rows to see how many dropped 
length(barcode9_dups) # 238 barcodes occur more than once
length(barcode9_unique) # 21,001 barcodes occur once
length(barcode9) - length(barcode9_unique) # 479 cells will be excluded

#export as tsv file since souporcell requires that
write_tsv((as.data.frame(barcode9_unique)), file = "~/pt9.barcodes.tsv", col_names = FALSE)

# Overlap between all combinations
P9_barcodes_ls <- P9 %>% select(barcode, orig.ident) %>%
  mutate(barcode = gsub("-.*", "", barcode)) %>%
  select(orig.ident, barcode) %>%
  group_by(orig.ident) %>% 
  pivot_wider(names_from = orig.ident, values_from = barcode, values_fn = list)

# Make a list of every sample compared to every other sample
overlap_list <- lapply(P9_barcodes_ls, function(x) lapply(P9_barcodes_ls, function(y) sum(duplicated(c(x[[1]], y[[1]])))))
overlap_list

# Generate a dataframe
do.call(cbind, overlap_list)
