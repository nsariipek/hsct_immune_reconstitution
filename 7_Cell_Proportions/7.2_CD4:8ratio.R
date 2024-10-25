
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

# Load the data that contains all cells+ assignments
seu <- readRDS(paste0(my_wd,"RDS files/assignment_seu.rds"))

#select only
meta= seu@meta.data  %>%
  dplyr::select(celltype,patient_identity,cohort, Sample,id, assignment) 
#choose donor or recipient cells if needed
  meta <- subset(x=meta, subset = assignment== "host")

meta <- subset(x=meta, subset = Sample %in% c("P01_1Rem", "P01_2Rem", "P02_1Rem", "P04_1Rem",  "P05_1Rem", "P06_1Rem", "P07_1Rem", "P08_1Rem"))
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
  mutate(cohort = gsub("cohort1", "Non-relapsed", gsub("cohort2", "Relapsed", cohort)),
cohort = factor(cohort, levels = c("Non-relapsed", "Relapsed")))%>%
  ggplot(aes(x=cohort, y=ratio)) +
  geom_jitter(aes(color=cohort), size=4)+
  ylab("CD4/CD8 Ratio") +
  theme_pubr()+
scale_color_manual(values = c("Relapsed"= "tomato1", "Non-relapsed"="royalblue1"))+
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
                     label.y= 4, tip.length = 1, size = 6, inherit.aes = TRUE) 
  

p1
pdf("CD4__CD8ratio.pdf", width = 5, height = 7.5)
p1
dev.off()
getwd()
