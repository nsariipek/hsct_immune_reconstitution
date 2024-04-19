# Nurefsan tries get some genes values across the patients in the cohorts 1-2

# Load the libraries
library(tidyverse)
library(Seurat)
library(ggforce)
library(ggplot2)
library(cowplot)
library(randomcoloR)
library(readxl)
library(data.table)
library(janitor)

# Clean the enviroment
rm(list=ls())

# For Nurefsan:
my_wd <- "/Users/dz855/Dropbox (Partners HealthCare)/ImmuneEscapeTP53/"

seu <- readRDS(paste0(my_wd, "AnalysisNurefsan/RDS files/seu_diet_merged.rds"))

# Add gene expression as metadata
metadata <- as_tibble(seu@meta.data, rownames = "cell")

mygene <- "ADGRG1"
# Add expression data to metadata
metadata$mygene <- LayerData(seu, layer = "data")[mygene,]

# Subset only 3-6 months remission samples
metadata.filter <- subset(x=metadata, subset = Sample %in% c("P01_1Rem", "P01_2Rem", "P02_1Rem", "P04_1Rem", "P05_1Rem", "P06_1Rem", "P07_1Rem", "P08_1Rem"))
# Compare cell number before and after filtering (check that it makes sense)
# Filter out samples with <10 cells
metadata.filter <- metadata.filter %>% group_by(Sample) %>% filter(n() >= 10)

# Barplot
p1 <- metadata.filter %>% group_by(Sample) %>%
  summarize(n = n(), mean_mygene = mean(mygene)) %>%
  ggplot(aes(x = Sample, y = mean_mygene)) +
  ylab("Mean normalized expression log(TP10K+1)") +
  geom_bar(stat="identity") +
  scale_x_discrete(drop=F) +
  ggtitle(mygene) +
  theme_bw() +
  theme(aspect.ratio = 0.5, axis.text.x = element_text(angle = 45, vjust= 1, hjust = 1, size = 10, color = "black"),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 12, color = "black"),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        plot.title = element_text(size = 14, hjust = 0.5))
p1

# Sina plot
p2 <- metadata.filter %>%
  ggplot(aes(x = Sample, y = mygene, color= cohort)) +
stat_sina() +
  geom_violin(scale = "width")+
  ylab("Normalized expression log(TP10K+1)") +
  ggtitle(mygene) +
  scale_x_discrete(drop = F) +
  theme_bw() +
  theme(aspect.ratio = 0.5,
        axis.text.x = element_text(angle = 45, vjust= 1, hjust = 1, size = 8, color = "black"),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 12, color = "black"),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        plot.title = element_text(size = 14, hjust = 0.5))

p2



