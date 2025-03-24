## Calculating CD4/8 Ratio
# NS 240418, updated at 250129

# Load libraries
library(tidyverse)
library(Seurat)
library(ggpubr)

# Empty environment
rm(list=ls())
# Set directory

setwd("/home/rstudio/TP53_ImmuneEscape/7_Cell_Proportions")

# Load the data that contains T cells+ assignments
Tcells <- readRDS("~/250128_Tcell_subset.rds")

# Select only needed variables
meta= Tcells@meta.data  %>%
      dplyr::select(celltype, cohort, sample_status, orig.ident, sample_id, patient_id,timepoint,survival, library_type) 
# Choose donor or recipient cells if needed
  #meta <- subset(x=meta, subset = assignment== "host")

rownames(meta) <- NULL

meta$type <- case_when(grepl("CD8.", meta$celltype)| grepl("γδ.", meta$celltype) ~ "CD8",
                       grepl("CD4.",meta$celltype)| grepl("Treg", meta$celltype)~ "CD4")
meta$type <- as.factor(meta$type)

#Add TP53 status
mt_patients <- paste0("P", sprintf("%02d", c(1:9, 10:12, 14, 17)))  # P01-P09 format
wt_patients <- paste0("P", sprintf("%02d", c(13, 15, 16, 18:33)))  # P13-P33 format

meta <- meta %>%
  mutate(patient_id = as.factor(patient_id))%>%
  mutate(TP53 = case_when(
    patient_id %in% mt_patients ~ "MT",
    patient_id %in% wt_patients ~ "WT"))

View(meta)

# Select only 100-day samples that were in remission at that timepoint
meta_subset<- meta %>%
             subset(sample_status == "remission" & timepoint %in% c("3", "5", "6") 
                    & TP53== "WT")

#For calculations make the table
tb <- 
  meta_subset %>% group_by(sample_id, type,cohort, survival,TP53) %>%
  dplyr::summarize(n = n())  %>%
  ungroup() %>%
  group_by(sample_id) %>% 
  pivot_wider(names_from = type, values_from = n) %>% mutate(ratio= CD4/CD8) 

  #View(tb)
  
#Visualize 
p1 <- tb %>%
  ggplot(aes(x=survival, y=ratio)) +
  geom_jitter(aes(color=TP53), size=4)+
  ylab("CD4/CD8 Ratio") +
  theme_pubr()+
 # scale_color_manual(values = c("Relapsed"= "tomato1", "Non-relapsed"="royalblue1"))+
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
  stat_compare_means(aes(group = survival), method = "t.test", 
                     method.args = list(var.equal = T),
                     label = "p.format", label.x = 1.6, 
                     label.y= 4, tip.length = 1, size = 6, inherit.aes = TRUE) 
  
p1
pdf("CD4__CD8ratio_percohort.pdf", width = 5, height = 7.5)
p1
dev.off()
getwd()
