## Calculating CD4/8 Ratio
# NS 250129 , updated at 240428 
# Load libraries
library(tidyverse)
library(Seurat)
library(ggpubr)

# Empty environment
rm(list=ls())

# Set directory
setwd("/home/rstudio/TP53_ImmuneEscape/5_Cell_Proportions")

# Load the seurat object and select and save the only T cells
seu <- readRDS("~/250426_Seurat_annotated.rds")

Tcells <-seu %>% subset(celltype %in% c("CD4 Naive", "CD4 Central Memory", "CD4 Effector Memory", "CD4 Regulatory", "CD8 Naive", "CD8 Central Memory", "CD8 Effector Memory 1", "CD8 Effector Memory 2", "CD8 Tissue Resident Memory", "T Proliferating"))

# Save this for future use   
saveRDS(Tcells,"~/250428_Tcells.rds")
       
# Levels disapear after turning seu object to metadata, add the new levels
my_levels <- c("CD4 Naive", "CD4 Central Memory", "CD4 Effector Memory", "CD4 Regulatory", "CD8 Naive", "CD8 Central Memory", "CD8 Effector Memory 1", "CD8 Effector Memory 2", "CD8 Tissue Resident Memory", "T Proliferating")    


# Select only needed variables
meta= Tcells@meta.data  %>%
      dplyr::select(celltype, cohort, sample_status, orig.ident, sample_id, patient_id,timepoint,cohort, library_type, TP53_status) 
# Choose donor or recipient cells if needed# Choose donor or recipient cells if needed cohort
  #meta <- subset(x=meta, subset = assignment== "host")

rownames(meta) <- NULL

meta$type <- case_when(grepl("CD8.", meta$celltype)| grepl("Gamma.", meta$celltype) ~ "CD8",
                       grepl("CD4.",meta$celltype)| grepl("Treg", meta$celltype)~ "CD4")
meta$type <- as.factor(meta$type)

# Select only 100-day samples that were in remission at that timepoint
meta_subset<- meta %>%
             subset(sample_status == "remission" & timepoint %in% c("3", "5", "6") 
              #  & TP53_status== "MUT"
              )

#For calculations make the table
tb <- 
  meta_subset %>% group_by(patient_id, type, cohort) %>%
  dplyr::summarize(n = n())  %>%
  ungroup() %>%
  group_by(patient_id) %>% 
  pivot_wider(names_from = type, values_from = n) %>% mutate(ratio= CD4/CD8) 

  #View(tb)

# Cohort colors
cohort_colors <- c("long-term-remission" = "#546fb5FF","relapse" = "#e54c35ff")
p1 <- tb %>%
  ggplot(aes(x = cohort, y = ratio)) +
  geom_jitter(aes(fill = cohort), shape = 21, size = 5, stroke = 0.2, width = 0.15, color = "black") +
  stat_summary(fun.data = mean_se, geom = "errorbar", 
               width = 0.3, linewidth = 0.6, color = "black") +
  stat_compare_means(
    aes(group = cohort),
    method = "t.test",
    method.args = list(var.equal = TRUE),
    label = "p.format",
    label.y = 3,
    size = 4,
    tip.length = 0.02) +
  scale_fill_manual(values = cohort_colors) +
  ylab("CD4/CD8 Ratio") +
  expand_limits(y = 0) +
  coord_cartesian(ylim = c(0, 3.5))+
  theme_minimal(base_size = 8) +
  theme(
    aspect.ratio = 2,
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 12),
    legend.position = "none",
    panel.grid = element_blank(),
    axis.ticks = element_line(color = "black"),  
    axis.ticks.length = unit(0.25, "cm"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4))

p1
# Save the plot
pdf("5.3_CD4-CD8ratio_all.pdf", width = 3, height = 6)
p1
dev.off()

