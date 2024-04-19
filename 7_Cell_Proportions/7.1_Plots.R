# Creating proportion graphs for each cohort
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
seu_diet_merged <- readRDS(paste0(my_wd,"AnalysisNurefsan/RDS files/seu_diet_merged.rds"))

# Alluvial plot

# Check the metadata
View(seu_diet_merged@meta.data) 

# For cohort 1  
t = seu_diet_merged@meta.data %>% rownames_to_column("barcode")%>% 
  filter(library_type=="MNC") %>% 
  filter(cohort == "cohort1")%>% 
  dplyr::select(groups,id,celltype,patient_identity,cohort) 

# See the new metadata 
head(t)

# Select the T cells
t <- subset(t, subset = celltype %in% c("CD4 Naïve","CD4 Memory","Treg","CD8 Naïve","CD8 Memory","CD8 Effector",
                                        "CD8 Terminally Exhausted","γδ T lymphocytes","NK T cells","CD56 Dim NK cells",
                                        "CD56 Bright NK cells"))

# Summarise and calculate the frequencies 
t_sum = t %>% 
  group_by(groups, cohort, celltype) %>%
  summarize(n = n()) %>%
  group_by(groups, cohort) %>%
  mutate(freq = n/sum(n)) %>%
  dplyr::select(-n)

head(t_sum)

#t cell order
ordered_celltypes = c("CD4 Naïve","CD4 Memory","Treg","CD8 Naïve","CD8 Memory","CD8 Effector",
                      "CD8 Terminally Exhausted","γδ T lymphocytes","NK T cells","CD56 Dim NK cells",
                      "CD56 Bright NK cells")
# all celltypes order
ordered_celltypes = c("CD4 Naïve","CD4 Memory","Treg","CD8 Naïve","CD8 Memory","CD8 Effector","CD8 Terminally Exhausted","γδ T lymphocytes","NK T cells","CD56 Dim NK cells","CD56 Bright NK cells","Pre B cells","B cells","Pro B cells","Plasma Cells",  "Early Erythroids","Mid Erythroids","Late Erythroids","Blasts","Monocytes","Non Classical Monocytes","Pro Monocytes","cDC","pDC","HSPCs","Unidentified","Doublets")

# re-order the groups
ordered_timepoints= c("pre-transplant", "remission(3-6mo)","remission(>6mo)","relapse")

# turn into a factor
t_sum$celltype = factor(t_sum$celltype, levels = ordered_celltypes)
t_sum$groups= factor(t_sum$groups, levels = ordered_timepoints)

# Plot 
p1 = ggplot(t_sum, aes(x =  as.factor(groups), fill = celltype, group = celltype,
                          stratum = celltype, alluvium = celltype, 
                          y = freq, label = celltype)) +  
   geom_stratum() + 
   geom_flow(stat = "alluvium") + 
   #geom_text(size=3, stat = "stratum", aes(label = after_stat(stratum))) + 
   theme(axis.text.x = element_text(face="plain", size=16, color="black"), 
        axis.title.x = element_text(face="plain", size=16, color="black"),
        axis.text.y = element_text(face="plain", size=12, color="black"),
        axis.title.y = element_blank(),
        plot.title = element_text(size=20, face="plain"),
       #legend.position = "none") +
       legend.text=element_text(size=18),legend.title=element_text(size=18)) +
  labs(x="Timepoints")
       
p1

# For cohort 2  
t = seu_diet_merged@meta.data %>% rownames_to_column("barcode")%>% 
  filter(library_type=="MNC") %>% 
  filter(cohort == "cohort2")%>% 
  dplyr::select(groups,id,celltype,patient_identity,cohort) 

# See the new metadata 
head(t)

# Select the T cells
t <- subset(t, subset = celltype %in% c("CD4 Naïve","CD4 Memory","Treg","CD8 Naïve","CD8 Memory","CD8 Effector",
                                        "CD8 Terminally Exhausted","γδ T lymphocytes","NK T cells","CD56 Dim NK cells",
                                        "CD56 Bright NK cells"))

# Summarise and calculate the frequencies 
t_sum = t %>% 
  group_by(groups, cohort, celltype) %>%
  summarize(n = n()) %>%
  group_by(groups, cohort) %>%
  mutate(freq = n/sum(n)) %>%
  dplyr::select(-n)

head(t_sum)

#t cell order
ordered_celltypes = c("CD4 Naïve","CD4 Memory","Treg","CD8 Naïve","CD8 Memory","CD8 Effector",
                      "CD8 Terminally Exhausted","γδ T lymphocytes","NK T cells","CD56 Dim NK cells",
                      "CD56 Bright NK cells")
# all celltypes order
ordered_celltypes = c("CD4 Naïve","CD4 Memory","Treg","CD8 Naïve","CD8 Memory","CD8 Effector","CD8 Terminally Exhausted","γδ T lymphocytes","NK T cells","CD56 Dim NK cells","CD56 Bright NK cells","Pre B cells","B cells","Pro B cells","Plasma Cells",  "Early Erythroids","Mid Erythroids","Late Erythroids","Blasts","Monocytes","Non Classical Monocytes","Pro Monocytes","cDC","pDC","HSPCs","Unidentified","Doublets")

# re-order the groups
ordered_timepoints= c("pre-transplant", "remission(3-6mo)","remission(>6mo)","relapse")

# turn into a factor
t_sum$celltype = factor(t_sum$celltype, levels = ordered_celltypes)
t_sum$groups= factor(t_sum$groups, levels = ordered_timepoints)

# Plot 
p2 = ggplot(t_sum, aes(x =  as.factor(groups), fill = celltype, group = celltype,
                       stratum = celltype, alluvium = celltype, 
                       y = freq, label = celltype)) +  
  geom_stratum() + 
  geom_flow(stat = "alluvium") + 
  #geom_text(size=3, stat = "stratum", aes(label = after_stat(stratum))) + 
  theme(axis.text.x = element_text(face="plain", size=16, color="black"), 
        axis.title.x = element_text(face="plain", size=16, color="black"),
        axis.text.y = element_text(face="plain", size=12, color="black"),
        axis.title.y = element_blank(),
        plot.title = element_text(size=20, face="plain"),
        legend.text=element_text(size=18),legend.title=element_text(size=18)) +
  labs(x="Timepoints")

p2


