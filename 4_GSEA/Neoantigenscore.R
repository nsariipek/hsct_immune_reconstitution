## Neantigen score+ TCRs

# Load the libraries
library(scRepertoire)
library(Seurat)
library(ggpubr)
library(tidyverse)
library(janitor)
library(rstatix)
library(Hmisc)
library(RColorBrewer)

# Empty environment
rm(list=ls())

# For Nurefsan:
my_wd <- "/Users/dz855/Dropbox (Partners HealthCare)/ImmuneEscapeTP53/"

# Load the seurat object from 6.1 script end of line 192 which has the TCR+ scRNA combined object 
combined <- readRDS(paste0(my_wd, "AnalysisNurefsan/RDS files/combined.RDS"))

# Add the neoantigen scores

# Add scores as a data table 
cd4neoA <- read.csv(paste0(my_wd,"AnalysisNurefsan/DGE/signatures/cd4.csv"), sep = "\t", header = T)
cd8neoA <- read.csv(paste0(my_wd,"AnalysisNurefsan/DGE/signatures/cd8.csv"), sep = "\t", header = T)
neoA <- read.csv(paste0(my_wd,"AnalysisNurefsan/DGE/signatures/neoantigen.csv"), sep = "\t", header = T)

#Subset for different celltypes
cd8cells <- subset(x = combined, subset = celltype %in% c("CD8 Effector","CD8 Memory","CD8 Naïve","CD8 Terminally Exhausted"))

cd8ef <- subset(x = combined, subset = celltype=="CD8 Effector")

# Add to the seurat object
neoantigen<- AddModuleScore(object = cd8ef, 
                            features = cd8neoA,
                            name = "neoantigen",
                            assay = "RNA",
                            search = T)

# View 
View(neoantigen@meta.data)

#Select only CD8+ T cells 
cd8 <- subset(x = neoantigen, subset = celltype %in% c("CD8 Effector","CD8 Memory","CD8 Naïve","CD8 Terminally Exhausted"))

cd8ef <- subset(x = neoantigen, subset = celltype=="CD8 Effector")
cd8mem <- subset(x = neoantigen, subset = celltype =="CD8 Memory")
cd8na <- subset(x = neoantigen, subset = celltype=="CD8 Naïve")
cd8ex <- subset(x = neoantigen, subset = celltype=="CD8 Terminally Exhausted")

cd4 <- subset(x = neoantigen, subset = celltype %in% c("CD4 Memory","CD4 Naïve","Treg"))

# Insanity check-Visualize it 
neoantigen <- SetIdent(neoantigen, value = "Sample")
VlnPlot(neoantigen, features = "neoantigen1", split.by = "cohort", sort = "increasing")

# turn into tibble
neotb <- as_tibble(neoantigen@meta.data, rownames = "cell") %>%
  select(Sample, CTnt, neoantigen1)

# For each TCR, calculate relative clonotype size (grouped by sample)
t <- neotb %>%
group_by(Sample) %>%
  mutate(n_total = n()) %>%
  ungroup() %>%
  group_by(Sample, CTnt) %>%
  dplyr::summarize(
    n = n(),
    prop = n / first(n_total),
    meanScore = mean(neoantigen1)) %>%
  ungroup()


p <- ggplot(t, aes(Sample, y=meanScore)) +
  geom_jitter(aes(size= prop))+
  theme(aspect.ratio = 1, axis.ticks.x = element_blank(), axis.text.x = element_blank())+
  theme_bw()


p
