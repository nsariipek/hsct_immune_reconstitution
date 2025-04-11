# Neantigen score+TCRs
# Nurefsan Sariipek, 24-07-01, updated at Terra 25-01-29

# Load the libraries
library(scRepertoire)
library(Seurat)
library(ggpubr)
library(tidyverse)
library(janitor)
library(rstatix)
library(RColorBrewer)
library(ggforce) # for geom_sina

# Empty environment
rm(list=ls())

#set wd
setwd("~/TP53_ImmuneEscape/3_DGE/3.1_Neoantigenscore/")

# Load the seurat object from 6.1 script end of line 163 which has the TCR+ scRNA combined object
combined <- readRDS("~/Tcells_TCR.rds")

# Define MT and WT patient groups using direct string matching
mt_patients <- paste0("P", sprintf("%02d", c(1:9, 10:12, 14, 17)))  # P01-P09 format
wt_patients <- paste0("P", sprintf("%02d", c(13, 15, 16, 18:33)))  # P13-P33 format

# Add the TP53 info to the data
combined@meta.data <- combined@meta.data %>%
  mutate(patient_id = as.factor(patient_id))%>%
  mutate(TP53 = case_when(
    patient_id %in% mt_patients ~ "MT",
    patient_id %in% wt_patients ~ "WT"))

# Load antigen scores from Rosenberg lab papers https://www.sciencedirect.com/science/article/pii/S1535610823003963?via%3Dihub
# https://www.science.org/doi/10.1126/science.abl5447?url_ver=Z39.88-2003&rfr_id=ori:rid:crossref.org&rfr_dat=cr_pub%20%200pubmed
cd4neoA <- read.csv("signatures/cd4.csv")
cd8neoA <- read.csv("signatures/cd8.csv")
neoA <- read.csv("signatures/neoantigen.csv")

# Optional
# #Only select the MT samples
# combined_subset_MT <- subset(x = combined, subset = TP53=="MT")
# #Only select the WT samples
# combined_subset_WT <- subset(x = combined, subset =TP53=="WT")

# Optional
# Subset for different celltypes before adding the antigen score
 cd8cells <- subset(x = combined, subset = celltype %in% c("CD8 Effector","CD8 Memory","CD8 Naïve","CD8 Exhausted"))
#  cd8ef <- subset(x = combined, subset = celltype == "CD8 Effector")
#  cd8mem <- subset(x = combined, subset = celltype =="CD8 Memory")
#  cd8na <- subset(x = combined, subset = celltype=="CD8 Naïve")
#  cd8ex <- subset(x = combined, subset = celltype=="CD8 Exhausted")
#  cd4cells <- subset(x = combined, subset = celltype %in% c("CD4 Memory","CD4 Naïve","Treg"))
#  treg <- subset(x = combined, subset = celltype == "Treg")
#  cd4mem <- subset(x = combined, subset = celltype == "CD4 Memory")
#  cd4na <- subset(x = combined, subset = celltype == "CD4 Naïve")

# Add this module score to the subsetted dataset
neoantigen <- AddModuleScore(object =   cd8cells,
                             features = cd8neoA,
                             name = "neoantigen",
                             assay = "RNA",
                             search = T)
# View 
#View(neoantigen@meta.data)
#you can see there are a lot of NAs

# If you want to select the t cell type after adding the antigen score Select only CD8+ T cells 
# neoantigen <- subset(x = neoantigen, subset = celltype %in% c("CD8 Effector","CD8 Memory","CD8 Naïve","CD8 Terminally Exhausted"))'
# neoantigen <- subset(x = neoantigen, subset = celltype == "CD8 Effector")
# neoantigen <- subset(x = neoantigen, subset = celltype =="CD8 Memory")
# neoantigen <- subset(x = neoantigen, subset = celltype=="CD8 Naïve")
# neoantigen <- subset(x = neoantigen, subset = celltype=="CD8 Terminally Exhausted")
# neoantigen <- subset(x = neoantigen, subset = celltype %in% c("CD4 Memory","CD4 Naïve","Treg"))

# Visualize the score per patient
neoantigen <- SetIdent(neoantigen, value = "sample_id")
p1 <- VlnPlot(neoantigen, features = "neoantigen1", split.by = "TP53", sort = "increasing")
p1

pdf("Neocd8cells_perpt_all.pdf", width = 20, height = 10)
p1
dev.off()

# Turn seurat object to tibble
neoantigen <- as_tibble(neoantigen@meta.data, rownames = "cell") 

# Select only needed variables
neotb <- neoantigen%>%
  select(sample_id, CTstrict, neoantigen1,survival,patient_id)

# For each TCR, calculate relative clonotype size (grouped by sample)
neotb_grouped <- neotb %>%
  group_by(sample_id) %>%
  mutate(n_total = n()) %>%
  ungroup() %>%
  group_by(sample_id, CTstrict,survival,patient_id) %>%
  dplyr::summarize(
    n = n(),
    prop = n / first(n_total),
    meanScore = mean(neoantigen1)) %>%
  ungroup()


# There may a difference between the cohorts when taking all clonotypes, but this test is not stringent enough
rem_meanScores <- filter(neotb_grouped, survival == "Non-relapsed")$meanScore
rel_meanScores <- filter(neotb_grouped, survival == "Relapsed")$meanScore
# T-test
t_test_result <- t.test(rem_meanScores, rel_meanScores)
# Wilcox test (more commonly used in the literature)
w_test_result <- wilcox.test(rem_meanScores, rel_meanScores)
# Fold change
fc <- median(rem_meanScores)/median(rel_meanScores)

# Sina plot, grouped by survival
s <- ggplot(neotb_grouped, aes(x=survival, y=meanScore)) +
  geom_sina(aes(size = prop, color = sample_id, group = survival), scale = "width")+
  geom_violin(alpha=0, scale = "width", draw_quantiles = 0.5) +
  #ggtitle(label = unique(neoantigen$celltype)) + # this should be changed in the final script
  theme_pubr() +
  theme(aspect.ratio = 1)+
  annotate("text", x = 1.5, y = max(neotb_grouped$meanScore)-0.02,
           label = paste0("Fold change: ", round(fc, 4), "\n",
             "p = ", signif(w_test_result$p.value, digits = 4)), size=4.5)

s

pdf("Neocd8cells_MTvsWT.pdf", width = 20, height = 10)
s
dev.off()

# Bar plot showing mean for each patient's score
b <-neotb_grouped %>%
  group_by(patient_id,survival) %>%
  summarize(mean_of_meanScore = mean(meanScore, na.rm = TRUE))%>%
  ggplot(aes(x = survival, y = mean_of_meanScore)) + 
  geom_bar(stat="summary", fun=mean, aes(fill= survival))+
  scale_fill_manual(values=survival_colors)+
  # scale_color_manual(values = c("orange", "#008080"))+
  geom_jitter(aes(color=patient_id), size = 4) +
  stat_summary(fun.data=mean_se, geom="errorbar", width=.5, linewidth=1) +
  theme_pubr() +
  coord_cartesian(ylim = c(0,1)) +
  theme(strip.text = element_text(size = 20 , color = "black", face="bold")) +
  theme(aspect.ratio = 2, axis.text.x = element_text(angle = 45, vjust= 1, hjust = 1,    size = 24, color = "black"),
        axis.title.x = element_blank(), 
        axis.text.y = element_text(size = 24),
        axis.title.y = element_text(size = 24, color = "black"),
        plot.title =  element_text(size=24,color="black", face="bold",hjust = 1 ),
        legend.key.size = unit(10,"mm"),
        legend.title = element_text(size = 20),
        legend.position = "right",
        legend.text = element_text(size = 20))

b




