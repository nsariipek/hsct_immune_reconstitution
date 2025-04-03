## Calculating CD4/8 Ratio
# NS 240418, updated at 250129

# Load libraries
library(tidyverse)
library(Seurat)
library(ggpubr)

# Empty environment
rm(list=ls())

# Set directory
setwd("/home/rstudio/TP53_ImmuneEscape/5_Cell_Proportions")

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
mt_patients <- paste0("P", sprintf("%02d", c(1:9, 10:12, 14, 17)))
wt_patients <- paste0("P", sprintf("%02d", c(13, 15, 16, 18:33)))  

meta <- meta %>%
  mutate(patient_id = as.factor(patient_id))%>%
  mutate(TP53 = case_when(
    patient_id %in% mt_patients ~ "MT",
    patient_id %in% wt_patients ~ "WT"))

View(meta)

# Select only 100-day samples that were in remission at that timepoint
meta_subset<- meta %>%
             subset(sample_status == "remission" & timepoint %in% c("3", "5", "6") 
                   & TP53== "MT"
                    )

#For calculations make the table
tb <- 
  meta_subset %>% group_by(sample_id, type,cohort, survival,TP53,patient_id) %>%
  dplyr::summarize(n = n())  %>%
  ungroup() %>%
  group_by(sample_id) %>% 
  pivot_wider(names_from = type, values_from = n) %>% mutate(ratio= CD4/CD8) 

  #View(tb)
  
p1 <- tb %>%
  mutate(survival = factor(survival, levels = c("Non-relapsed","Relapsed"))) %>%
  ggplot(aes(x = survival, y = ratio)) +
  geom_jitter(aes(fill = survival), shape = 21, size = 5, stroke = 0.2, width = 0.15, color = "black") +
  stat_summary(fun.data = mean_se, geom = "errorbar", 
               width = 0.3, linewidth = 0.6, color = "black") +
  stat_compare_means(
    aes(group = survival),
    method = "t.test",
    method.args = list(var.equal = TRUE),
    label = "p.signif",
    label.y = 3,
    size = 4,
    tip.length = 0.02) +
  scale_fill_manual(values = c("Relapsed" = "#E64B35FF", "Non-relapsed" = "#4DBBD5FF")) +
  ylab("CD4/CD8 Ratio") +
  expand_limits(y = 0) +
  coord_cartesian(ylim = c(0, 3.5))+
  theme_minimal(base_size = 8) +
  theme(
    aspect.ratio = 1,
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 12),
    legend.position = "none",
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4))

p1
pdf("5.3_CD4-CD8ratio_MT_only.pdf", width = 3, height = 3)
p1
dev.off()
getwd()
