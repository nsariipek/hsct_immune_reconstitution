# Peter van Galen, 250710

# Load the libraries
library(tidyverse)

# Empty environment
rm(list = ls())

# For Peter
setwd(
  "~/DropboxMGB/Projects/ImmuneEscapeTP53/hsct_immune_reconstitution/05_DGE"
)

# Load the DGE results from 5.5_DGE_tumorcells.R
de_results <- read_tsv("5.5_DGE_Pre-transplant_vs_Relapse.txt")

# MHC-II genes. This is probably not exhaustive.
genes <- c(
  "CIITA",
  "HLA-DMA",
  "HLA-DMB",
  "HLA-DOA",
  "HLA-DOB",
  "HLA-DPA1",
  "HLA-DPB1",
  "HLA-DQA1",
  "HLA-DQA2",
  "HLA-DQB1",
  "HLA-DQB1-AS1",
  "HLA-DQB2",
  "HLA-DRA",
  "HLA-DRB1",
  "HLA-DRB5",
  "CD74"
)
# None of these are significant

# Toffalori et al.
# "We confirmed the deregulation of multiple costimulatory ligands on AML blasts at post-transplantation relapse (PD-L1, B7-H3, CD80, PVRL2)"
genes <- c("CD274", "CD276", "CD80", "NECTIN2")
# Of these, only NECTIN2 is up in our data

de_results %>% filter(gene %in% genes)
