# Creating proportion graphs for each cohort as alluvial plots
# Updated by Nurefsan at 240419
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
library(gridExtra)
library(cowplot)
library(patchwork)

# Empty environment
rm(list=ls())

# For Nurefsan:
my_wd <- "/Users/dz855/Dropbox (Partners HealthCare)/ImmuneEscapeTP53/"
# For Peter:
#my_wd <- "~/DropboxMGB/Projects/ImmuneEscapeTP53/"

# Load the data
seu_diet_merged <- readRDS(paste0(my_wd,"/RDS files/seu_diet_merged.rds"))

# Check the metadata
#View(seu_diet_merged@meta.data) 

############ Check the M/L ratio at 100 day in cohorts 1-2 ###############################

tbl = seu_diet_merged@meta.data %>% rownames_to_column("barcode") %>% 
  filter(library_type=="MNC") %>% 
  filter(cohort %in% c("cohort1", "cohort2")) %>% 
  mutate(celltype= gsub("Blast","blast",celltype))
  
tbl_filtered <- subset(x = tbl, subset = Sample %in% c("P01_1Rem", "P01_2Rem", "P02_1Rem", "P04_1Rem", "P05_1Rem", "P06_1Rem", "P07_1Rem", "P08_1Rem"))


# Add a new column naming if it is a L or M cell, L = T, NK, B cells

tbl_filtered$type <- case_when(
  grepl("CD|γδ|NK|Plasma|B|Treg\\b", tbl_filtered$celltype) ~ "L",
  grepl("Monocytes|Erythroids|DC", tbl_filtered$celltype) ~ "M",
  TRUE ~ NA_character_  # Default case if none of the above conditions match
)                           
                             
# check to see if this work properly 
tbl_filtered %>% tabyl(celltype, type)

tbl2 <-
  tbl_filtered %>% group_by(Sample, type,cohort) %>%
  dplyr::summarize(n = n())  %>%
  ungroup() %>%
  group_by(Sample) %>% 
  pivot_wider(names_from = type, values_from = n) %>% mutate(ratio = M/L)

tbl2 %>% ggplot(aes(x = cohort, y = ratio)) +
  geom_jitter(width = 0.1)
  
#####################################################
# Summarise the data for cohort 1  
t = seu_diet_merged@meta.data %>% rownames_to_column("barcode")%>% 
  filter(library_type=="MNC") %>% 
  filter(cohort == "cohort1")%>% 
  dplyr::select(groups,id,celltype,patient_identity,cohort) 

# See the new metadata 
head(t)

# Select the T cells only
t <- subset(t, subset = celltype %in% c("CD4 Naïve","CD4 Memory","Treg","CD8 Naïve","CD8 Memory","CD8 Effector", "CD8 Terminally Exhausted","γδ T lymphocytes","NK T cells","CD56 Dim NK cells","CD56 Bright NK cells"))

# Summarise and calculate the frequencies 
t_sum = t %>% 
  group_by(groups, cohort, celltype) %>%
  summarize(n = n()) %>%
  group_by(groups, cohort) %>%
  mutate(freq = n/sum(n)) %>%
  dplyr::select(-n)

head(t_sum)

# Order the T cells from immature to mature
ordered_celltypes = c("CD4 Naïve","CD4 Memory","Treg","CD8 Naïve","CD8 Memory","CD8 Effector","CD8 Terminally Exhausted","γδ T lymphocytes","NK T cells","CD56 Bright NK cells","CD56 Dim NK cells")
# all celltypes order
# ordered_celltypes = c("CD4 Naïve","CD4 Memory","Treg","CD8 Naïve","CD8 Memory","CD8 Effector","CD8 Terminally Exhausted","γδ T lymphocytes","NK T cells","CD56 Dim NK cells","CD56 Bright NK cells","Pre B cells","B cells","Pro B cells","Plasma Cells",  "Early Erythroids","Mid Erythroids","Late Erythroids","Blasts","Monocytes","Non Classical Monocytes","Pro Monocytes","cDC","pDC","HSPCs","Unidentified","Doublets")

# re-order the groups to make sense on the graph
ordered_timepoints= c("pre-transplant", "remission(3-6mo)","remission(>6mo)","relapse")

# turn both vectors into a factor
t_sum$celltype = factor(t_sum$celltype, levels = ordered_celltypes)
t_sum$groups= factor(t_sum$groups, levels = ordered_timepoints)


# Add the color palette that Peter has used in his UMAP in code 6.1
qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pdf("230724_Colorwheel.pdf")
pie(rep(1,74), col=col_vector)
dev.off()

mycol_tib <- tribble(~celltype, ~color,
                     "CD4 Naïve", col_vector[23],
                     "CD4 Memory", col_vector[28],
                     "Treg", col_vector[11],
                     "CD8 Naïve", col_vector[67],
                     "CD8 Memory", col_vector[63],
                     "CD8 Effector", col_vector[58],
                     "CD8 Terminally Exhausted", col_vector[60],
                     "γδ T lymphocytes", col_vector[1],
                     "NK T cells", col_vector[46],
                     "CD56 Bright NK cells", col_vector[61],
                     "CD56 Dim NK cells", col_vector[52])
mycol <- mycol_tib$color
names(mycol) <- mycol_tib$celltype

# Plot 1, Cohort 1
p1 = ggplot(t_sum, aes(x =  as.factor(groups), fill = celltype, group = celltype,
                          stratum = celltype, alluvium = celltype, 
                          y = freq, label = celltype)) +  
  geom_stratum() + 
  geom_flow(stat = "alluvium") + 
  scale_fill_manual(values = mycol) + 
  theme_bw()+
   #geom_text(size=3, stat = "stratum", aes(label = after_stat(stratum))) + 
   theme(aspect.ratio = 0.75,
        axis.text.x = element_text(face="plain", size=16, color="black"), 
        axis.title.x = element_text(face="plain", size=16, color="black"),
        axis.title.y = element_text(face="plain", size=16, color="black"),
        axis.text.y = element_text(face="plain", size=12, color="black"),
        plot.title = element_text(size=20, face="plain"),
       #legend.position = "none") +
       legend.text=element_text(size=18),legend.title=element_text(size=18)) +
       labs(x="Timepoints", y="Proportions")
       
p1

# Summarise the data for cohort 2  
y = seu_diet_merged@meta.data %>% rownames_to_column("barcode")%>% 
  filter(library_type=="MNC") %>% 
  filter(cohort == "cohort2")%>% 
  dplyr::select(groups,id,celltype,patient_identity,cohort) 

# See the new metadata 
head(y)

# Select the T cells only
y <- subset(y, subset = celltype %in% c("CD4 Naïve","CD4 Memory","Treg","CD8 Naïve","CD8 Memory","CD8 Effector",
                                        "CD8 Terminally Exhausted","γδ T lymphocytes","NK T cells","CD56 Bright NK cells","CD56 Dim NK cells"
                                        ))
# Summarise and calculate the frequencies 
y_sum = y %>% 
 group_by(groups, cohort, celltype) %>%
  summarize(n = n()) %>%
  group_by(groups, cohort) %>%
  mutate(freq = n/sum(n)) %>%
  dplyr::select(-n)

head(y_sum)

# turn into a vectors to factor to re-order it with saved characters from Cohort1
y_sum$celltype = factor(y_sum$celltype, levels = ordered_celltypes)
y_sum$groups= factor(y_sum$groups, levels = ordered_timepoints)

# Plot2, Cohort 2
p2 = ggplot(y_sum, aes(x =  as.factor(groups), fill = celltype, group = celltype,
                       stratum = celltype, alluvium = celltype, 
                       y = freq, label = celltype)) +  
  geom_stratum() + 
  geom_flow(stat = "alluvium") + 
  scale_fill_manual(values = mycol) + 
  theme_bw()+
  #geom_text(size=3, stat = "stratum", aes(label = after_stat(stratum))) + 
  theme(aspect.ratio=0.75,
          axis.text.x = element_text(face="plain", size=16, color="black"), 
        axis.title.x = element_text(face="plain", size=16, color="black"),
        axis.text.y = element_text(face="plain", size=12, color="black"),
        axis.title.y = element_blank(),
        plot.title = element_text(size=20, face="plain"),
        legend.position = "none",
        legend.text=element_text(size=18),legend.title=element_text(size=18)) +
  labs(x="Timepoints")

p2


# Save the plots in one file
p <- plot_grid(p1,p2)
save_plot("cohortscombined.pdf", p, ncol = 2, base_asp = 2, base_height = 8, base_width = 16)

#or indivually
#p1
pdf("Cohort1.pdf", width = 16, height = 8)
p1
dev.off()

#p2
pdf("Cohort2.pdf", width = 16, height = 8)
p2
dev.off()

# to see
#p1+p2


