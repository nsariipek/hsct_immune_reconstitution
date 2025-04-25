# Neoantigen score+TCRs
# Nurefsan Sariipek, 24-07-01, updated at 25-04-21
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

# Set working directory
setwd("~/TP53_ImmuneEscape/3_DGE/3.3_Neoantigenscore/")

# Load the seurat object from 6.1 script end of line 163 which has the TCR+ scRNA combined object
combined <- readRDS("~/250421_Tcells_TCR.rds")

# Select only 3-6 mo and remission samples
combined_subset <- subset(x= combined, subset= timepoint %in% c("3","5","6") & sample_status == "remission")

# Load antigen scores from Rosenberg lab papers https://www.sciencedirect.com/science/article/pii/S1535610823003963?via%3Dihub
# https://www.science.org/doi/10.1126/science.abl5447?url_ver=Z39.88-2003&rfr_id=ori:rid:crossref.org&rfr_dat=cr_pub%20%200pubmed
cd4neoA <- read.csv("signatures/cd4.csv")
cd8neoA <- read.csv("signatures/cd8.csv")
neoA <- read.csv("signatures/neoantigen.csv")

# Optional step for subsetting
# Subset for different celltypes before adding the antigen score
 cd8cells <- subset(x = combined_subset, subset = celltype %in% c("CD8 Effector","CD8 Memory","CD8 Naive","CD8 Exhausted","Delta-Gamma T"))
#  cd8ef <- subset(x = combined_subset, subset = celltype == "CD8 Effector")
#  cd8mem <- subset(x = combined_subset, subset = celltype =="CD8 Memory")
#  cd8na <- subset(x = combined_subset, subset = celltype=="CD8 Naïve")
#  cd8ex <- subset(x = combined, subset = celltype=="CD8 Exhausted")
#  cd4cells <- subset(x = combined, subset = celltype %in% c("CD4 Memory","CD4 Naïve","Treg"))
#  treg <- subset(x = combined, subset = celltype == "Treg")
#  cd4mem <- subset(x = combined, subset = celltype == "CD4 Memory")
#  cd4na <- subset(x = combined, subset = celltype == "CD4 Naïve")

# Add this module score to the subsetted dataset
neoantigen <- AddModuleScore(object = combined_subset,
                             features = neoA,
                             name = "neoantigen",
                             assay = "RNA",
                             search = T)

# Visualize the score per patient
neoantigen <- SetIdent(neoantigen, value = "patient_id")
p1 <- VlnPlot(neoantigen, features = "neoantigen1", split.by = "TP53", sort = "increasing")
p1

pdf("Neocd8cells_perpt_all.pdf", width = 20, height = 10)
p1
dev.off()

# Turn seurat object to tibble
neoantigen <- as_tibble(neoantigen@meta.data, rownames = "cell") 

# Select only needed variables
neotb <- neoantigen %>%
  select(patient_id,CTstrict, neoantigen1,cohort,patient_id,TP53_status)

# For each TCR, calculate relative clonotype size (grouped by sample)
neotb_grouped <- neotb %>% 
  filter(cohort=="long-term-remission") %>%
  group_by(patient_id) %>%
  mutate(n_total = n()) %>%
  ungroup() %>%
  group_by (CTstrict,TP53_status,cohort,patient_id) %>%
  dplyr::summarize(
    n = n(),
    prop = n / first(n_total),
    meanScore = mean(neoantigen1)) %>%
    ungroup()

# There may a difference between the cohorts when taking all clonotypes, but this test is not stringent enough
mt_meanScores <- filter(neotb_grouped, TP53_status == "MUT")$meanScore
wt_meanScores <- filter(neotb_grouped, TP53_status == "WT")$meanScore
# T-test
t_test_result <- t.test(mt_meanScores, wt_meanScores)
# Wilcox test (more commonly used in the literature)
w_test_result <- wilcox.test(mt_meanScores, wt_meanScores)
# Fold change
fc <- median(mt_meanScores)/median(wt_meanScores)

# Sina plot, grouped by survival
s <- ggplot(neotb_grouped, aes(x=TP53_status, y=meanScore)) +
  geom_sina(aes(size = prop, color = patient_id, group = TP53_status), scale = "width")+
  geom_violin(alpha=0, scale = "width", draw_quantiles = 0.5) +
  #ggtitle(label = unique(neoantigen$celltype)) + # this should be changed in the final script
  theme_pubr() +
  theme(aspect.ratio = 1)+
  annotate("text", x = 1.5, y = max(neotb_grouped$meanScore)-0.02,
           label = paste0("Fold change: ", round(fc, 4), "\n",
             "p = ", signif(w_test_result$p.value, digits = 8)), size=4.5)
s


pdf("Relapse_cohort_Neo_cells_MTvsWT.pdf", width = 20, height = 10)
s
dev.off()

# survival_colors <- c("long-term-remission" = "#546fb5FF","relapse" = "#e54c35ff")
# # Bar plot showing mean for each patient's score
# b <-neotb_grouped %>%
#   group_by(patient_id,survival) %>%
#   summarize(mean_of_meanScore = mean(meanScore, na.rm = TRUE))%>%
#   ggplot(aes(x = survival, y = mean_of_meanScore)) + 
#   geom_bar(stat="summary", fun=mean, aes(fill= survival))+
#   scale_fill_manual(values=survival_colors)+
#   # scale_color_manual(values = c("orange", "#008080"))+
#   geom_jitter(aes(color=patient_id), size = 4) +
#   stat_summary(fun.data=mean_se, geom="errorbar", width=.5, linewidth=1) +
#   theme_pubr() +
#   coord_cartesian(ylim = c(0,1)) +
#   theme(strip.text = element_text(size = 20 , color = "black", face="bold")) +
#   theme(aspect.ratio = 2, axis.text.x = element_text(angle = 45, vjust= 1, hjust = 1,    size = 24, color = "black"),
#         axis.title.x = element_blank(), 
#         axis.text.y = element_text(size = 24),
#         axis.title.y = element_text(size = 24, color = "black"),
#         plot.title =  element_text(size=24,color="black", face="bold",hjust = 1 ),
#         legend.key.size = unit(10,"mm"),
#         legend.title = element_text(size = 20),
#         legend.position = "right",
#         legend.text = element_text(size = 20))
# 
# b




