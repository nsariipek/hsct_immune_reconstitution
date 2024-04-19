
## Calculating CD4/8 Ratio
#Updated by NS 240418

# Load libraries
library(tidyverse)
library(Seurat)
library(ggplot2)
library(randomcoloR)
library(readxl)
library(data.table)
library(ggforce)
library("RColorBrewer")
library(ggpubr)
library(janitor)
library(purrr)
library(ggrepel)
library(ggalluvial)

# Empty environment
rm(list=ls())

# For Nurefsan:
my_wd <- "/Users/dz855/Dropbox (Partners HealthCare)/ImmuneEscapeTP53/"

# Load the data
Tcells <- readRDS(paste0(my_wd,"AnalysisNurefsan/RDS files/Tcellsfinal.rds"))

#select only
meta= Tcells@meta.data  %>%
  dplyr::select(groups,celltype,patient_identity,cohort, Sample,id)


meta <- subset(x=meta, subset = id %in% c("P01_1Rem", "P01_1RemT", "P01_2Rem", "P02_1Rem", "P04_1Rem", "P04_1RemT", "P05_1Rem", "P06_1Rem", "P07_1Rem", "P07_1RemT", "P08_1Rem", "P08_1RemT"))
rownames(meta) <- NULL

meta$type <- case_when(grepl("CD8.", meta$celltype) ~ "cd8",
                       grepl("CD4.",meta$celltype)| grepl("Treg", meta$celltype)~ "cd4",
                       grepl("γδ.",  meta$celltype)| grepl("NK", meta$celltype)  ~ " others")
meta$type <- as.factor(meta$type)
View(meta)

t2 <- 
  meta %>% group_by(Sample, type,cohort) %>%
 dplyr::summarize(n = n())  %>%
  ungroup() %>%
  group_by(Sample) %>% 
  pivot_wider(names_from = type, values_from = n) %>% mutate(ratio= cd4/cd8) 

  #rownames(t2) <- t2$type
  #t2$type <- NULL
  
View(t2)

p1 <- t2 %>%
  ggplot(aes(x=cohort, y=ratio)) +
  geom_jitter(aes(color=cohort), size=4)+
  theme_bw() +
  ylab("CD4/CD8 Ratio") +
  theme_pubr()+
  theme(strip.text = element_text(size = 14, color = "black", face="bold"),
        aspect.ratio = 1.5, 
        axis.text.x = element_text(angle = 45, vjust= 1, hjust = 1, size = 15, color = "black"),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 15, color = "black"),
        legend.key.size = unit(3,"mm"),
        legend.position = "right",
        legend.text = element_text(size = 12)) +
  expand_limits(x = 0, y = 0) + 
  stat_summary(fun.data=mean_se, geom="errorbar", width=.5, linewidth=1) +
  stat_compare_means(aes(group = cohort), method = "t.test", 
                     method.args = list(var.equal = T),
                     label = "p.format", label.x = 1.6, 
                     label.y= 2, tip.length = 1, size = 6, inherit.aes = TRUE) 
  

p1
pdf("CD4-CD8ratiopercohort.pdf", width = 5, height = 7.5)
p1
dev.off()
getwd()
