# Peter van Galen, 241101
# Visualize one of the Numbat results files 
# This code has not been updated to work in this location

library(tidyverse)
library(igraph)
library(numbat)
library(glue)
library(data.table)
library(ggtree)
library(tidygraph)
library(patchwork)

# Empty environment
rm(list=ls())

# For Nurefsan:
setwd("~/TP53_ImmuneEscape/5_Cell_Proportions/")

# Favorite function
cutf <- function(x, f=1, d="/") sapply(strsplit(x, d), function(i) paste(i[f], collapse=d))

# Test: download example data via http://pklab.med.harvard.edu/teng/data/nb_TNBC1.rds
#nb <- readRDS("~/Library/Mobile Documents/com~apple~CloudDocs/Downloads/nb_TNBC1.rds")

# Load numbat object from our own data (need to mount Broad storage first)
nb = Numbat$new(out_dir = "/Volumes/broad_vangalenlab/sariipek/numbat/P05")

# Generate CNV plot
mypal = c('1' = 'gray', '2' = "#377EB8", '3' = "#4DAF4A", '4' = "#984EA3")

nb$plot_phylo_heatmap(
  clone_bar = TRUE, 
  p_min = 0.9,
  pal_clone = mypal
)


