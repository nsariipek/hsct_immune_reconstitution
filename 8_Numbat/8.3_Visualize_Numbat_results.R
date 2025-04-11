# Peter van Galen, 241101
# Updated by Nurefsan Sariipek, 250411
# Visualize one of the Numbat results files 

# Better to run this script on the local space rather than Terra
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
setwd("/Users/dz855/Partners HealthCare Dropbox/Nurefsan Sariipek/ImmuneEscapeTP53/TP53_ImmuneEscape/8_Numbat")

# Favorite function
cutf <- function(x, f=1, d="/") sapply(strsplit(x, d), function(i) paste(i[f], collapse=d))

# Test: download example data via http://pklab.med.harvard.edu/teng/data/nb_TNBC1.rds
#nb <- readRDS("~/Library/Mobile Documents/com~apple~CloudDocs/Downloads/nb_TNBC1.rds")

# Load numbat object from our own data (need to mount Broad storage first)
#For Nurefsan
nb = Numbat$new(out_dir = "/Volumes/sariipek/numbat/old/P08")
pagoda = readRDS(url('http://pklab.org/teng/data/con_TNBC1.rds'))


# Generate CNV plot
mypal = c('1' = 'gray', '2' = "#377EB8", '3' = "#4DAF4A", '4' = "#984EA3")

plot1 <- nb$plot_phylo_heatmap(
  clone_bar = TRUE, 
  p_min = 0.9,
  pal_clone = mypal
)
#Save the panel as pdf
plot1
pdf("P08_panel2.pdf", width=15, height=5)
plot1
dev.off()

plots = lapply(
  1:4,
  function(k) {
    nb$cutree(n_cut = k)
    nb$plot_phylo_heatmap() + ggtitle(paste0('n_cut=', k))
  }
)
wrap_plots(plots)

# Visualize  CNV events in pseudobulks where the data is more rich, aggregating cells by clone
plot2 <- nb$bulk_clones %>% 
  filter(n_cells > 50) %>%
  plot_bulks(
    min_LLR = 10, # filtering CNVs by evidence
    legend = TRUE
  )
plot2
#Save as pdf
pdf("P08_bulk_clones.pdf", width=15, height=10)
plot2
dev.off()



