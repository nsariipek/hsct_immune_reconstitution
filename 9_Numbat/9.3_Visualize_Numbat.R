# Peter van Galen, 241101
# Visualize one of the Numbat results files 
# This code has not been updated to work in this location

library(tidyverse)
library(igraph)
library(numbat)
library(numbat)
library(glue)
library(data.table)
library(ggtree)
library(tidygraph)
library(patchwork)

# Favorite function
cutf <- function(x, f=1, d="/") sapply(strsplit(x, d), function(i) paste(i[f], collapse=d))

# Load numbat object. First download example data via http://pklab.med.harvard.edu/teng/data/nb_TNBC1.rds
nb <- readRDS("nb_TNBC1.rds")

# Load numbat object from our own data (need to mount Broad storage first)
nb = Numbat$new(out_dir = '/Volumes/broad_vangalenlab/sariipek/numbat/2737/host_2737')

# Generate CNV plot
mypal = c('1' = 'gray', '2' = "#377EB8", '3' = "#4DAF4A", '4' = "#984EA3")

nb$plot_phylo_heatmap(
  clone_bar = TRUE, 
  p_min = 0.9,
  pal_clone = mypal
)


# BARCODE TROUBLESHOOTING (THIS IS OBSOLETE)

bcs <- read_tsv("~/Library/Mobile Documents/com~apple~CloudDocs/Downloads/Numbat_2737_pt5_rel_host_barcodes.tsv", col_names = "cell")
length(nb$clone_post$cell)

# All input barcodes are reported by Numbat Step 2
all( bcs$cell %in% nb$clone_post$cell )

# But some barcodes are "added" by Numbat Step 2...!
mystery_barcodes <- setdiff(nb$clone_post$cell, bcs$cell)

# Get souporcell calls
my_wd <- "~/DropboxMGB/Projects/ImmuneEscapeTP53/TP53_ImmuneEscape/5_Souporcell/"
combined_df <- read_csv(paste0(my_wd,file = "/results/cohort1-2_souporcell.csv"))
pt5 <- subset(combined_df, subset = patient_identity=="pt05")

# This validates the length of the barcode list provided to numbat
pt5 %>% filter(status == "relapse", assignment == "host", library_type == "MNC") # 8,247
length( bcs$cell )

# Why does Numbat Step 2 create more barcodes than it was provided with???
mystery_barcodes %in% bcs$cell
# 709 are from pt5
sum( mystery_barcodes %in% cutf(pt5$cell, d = "_") )
# 809 are in the combined object
sum( mystery_barcodes %in% cutf(combined_df$cell, d = "_") )

# What about the provided counts?
testRDS <- readRDS("/Volumes/broad_vangalenlab/sariipek/numbat/2737/pt5_rel_host_numbat_umi_counts.rds")
dim(testRDS)
# So we found the problem: the provided rds object has too many cells