# Peter and Nurefsan, 241025
# Figure out why the Numbat results looked very similar whether or not the cell barcode list was subsetted for host.

library(tidyverse)

cutf <- function(x, f=1, d="/") sapply(strsplit(x, d), function(i) paste(i[f], collapse=d))

# Load barcodes that were used to start Numbat
# First download them (zsh)
# cd ~
# gsutil cp gs://ns_bucket_220819/Numbat/2737/2737_barcodes.tsv.gz .
# gsutil cp gs://ns_bucket_220819/Numbat/2737/pt5_rel_host_barcodes.tsv .
# Load in R
all_bcs <- read_tsv("~/2737_barcodes.tsv", col_names = "BC")
host_bcs <- read_tsv("~/pt5_rel_host_barcodes.tsv", col_names = "BC")


# STEP 1 CHECKS

# Load Numbat Step 1 output
numbat_step1_all <- read_tsv("/Volumes/broad_vangalenlab/sariipek/numbat/2737/2737_allele_counts.tsv")
numbat_step1_host <- read_tsv("/Volumes/broad_vangalenlab/sariipek/numbat/2737/2737_host_allele_counts.tsv")

# We provided 9517 barcodes as input
length(all_bcs$BC)
# Step 1 created data for 9447 barcodes
length(unique(numbat_step1_all$cell))
# All of these are in the input barcodes
all( unique(numbat_step1_all$cell) %in% all_bcs$BC )

# We provided 7728 barcodes as input
length(host_bcs$BC)
# Step 1 created data for 7728 barcodes
length(unique(numbat_step1_host$cell))
# All of these are in the input barcodes
all( unique(numbat_step1_host$cell) %in% host_bcs$BC )

# This strongly suggests Step 1 filtered for cell barcodes as intended


# STEP 2 CHECKS

# Load Numbat Step 2 output
numbat_step2_all1 <- read_tsv("/Volumes/broad_vangalenlab/sariipek/numbat/2737/geno_1.tsv")
numbat_step2_all2 <- read_tsv("/Volumes/broad_vangalenlab/sariipek/numbat/2737/geno_2.tsv")
numbat_step2_host <- read_tsv("/Volumes/broad_vangalenlab/sariipek/numbat/2737/host_2737/geno_1.tsv")

# Surprisingly, the Step 2 output has the same 8941 cell barcodes in both analysis (all & host)
identical(numbat_step2_all1$cell, numbat_step2_all2$cell)
identical(cutf(numbat_step2_all1$cell, d = "_", f = 1), numbat_step2_host$cell)

# Compare with number of cells in the umi_counts file: also 8941
umiCounts <- readRDS("/Volumes/broad_vangalenlab/sariipek/numbat/2737/2737_numbat_umi_counts.rds")
dim(umiCounts)

# We conclude that Numbat Step 2 returns results based on the count matrix, so the count matrix should be subsetted for Numbat Step 2 if we only want to analyze host cells


# Check: how can numbat_step2_host have data for cells that were not in host_bcs? This is somewhat confusing.
numbat_step2_host %>% filter(! cell %in% host_bcs$BC)

