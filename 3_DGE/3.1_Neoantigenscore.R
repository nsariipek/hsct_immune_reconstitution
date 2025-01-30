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

# Load antigen scores from Rosenberg lab papers https://www.sciencedirect.com/science/article/pii/S1535610823003963?via%3Dihub, 
# https://www.science.org/doi/10.1126/science.abl5447?url_ver=Z39.88-2003&rfr_id=ori:rid:crossref.org&rfr_dat=cr_pub%20%200pubmed
cd4neoA <- read.csv("signatures/cd4.csv")
cd8neoA <- read.csv("signatures/cd8.csv")
neoA <- read.csv("signatures/neoantigen.csv")

# Subset for different celltypes before adding the antigen score
# cd8cells <- subset(x = combined, subset = celltype %in% c("CD8 Effector","CD8 Memory","CD8 Naïve","CD8 Exhausted"))
#  cd8ef <- subset(x = combined, subset = celltype == "CD8 Effector")
#  cd8mem <- subset(x = combined, subset = celltype =="CD8 Memory")
#  cd8na <- subset(x = combined, subset = celltype=="CD8 Naïve")
#  cd8ex <- subset(x = combined, subset = celltype=="CD8 Exhausted")
cd4cells <- subset(x = combined, subset = celltype %in% c("CD4 Memory","CD4 Naïve","Treg"))
treg <- subset(x = combined, subset = celltype == "Treg")
cd4mem <- subset(x = combined, subset = celltype == "CD4 Memory")
cd4na <- subset(x = combined, subset = celltype == "CD4 Naïve")

# Add to the seurat object
neoantigen <- AddModuleScore(object =   cd8cells,
                             features = neoA,
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

# Sanity check, visualize the violin plots
neoantigen <- SetIdent(neoantigen, value = "sample_id")
VlnPlot(neoantigen, features = "neoantigen1", split.by = "survival", sort = "increasing")

# Turn into tibble
neoantigen <- as_tibble(neoantigen@meta.data, rownames = "cell") 

neotb <- neoantigen%>%
  select(sample_id, CTstrict, neoantigen1,survival)

# For each TCR, calculate relative clonotype size (grouped by sample)
neotb_grouped <- neotb %>%
  group_by(sample_id) %>%
  mutate(n_total = n()) %>%
  ungroup() %>%
  group_by(sample_id, CTstrict,survival) %>%
  dplyr::summarize(
    n = n(),
    prop = n / first(n_total),
    meanScore = mean(neoantigen1)) %>%
  ungroup()
# More elegant version should yield the exact same result:
# neotb_grouped <- neotb %>%
#   group_by(Sample, CTstrict) %>%
#   dplyr::summarize(
#     n = n(),
#     meanScore = mean(neoantigen1)) %>%
#   ungroup() %>%
#   group_by(Sample) %>%
#   mutate(prop = n / sum(n))

# med_df <- aggregate(meanScore ~ Sample, data = neotb_grouped, FUN = median)
# 
# 
# p <- ggplot(neotb_grouped, aes(x=Sample, y=meanScore)) +
#   geom_jitter(aes(size= prop), width = 0.25)+
#   geom_crossbar(data= med_df, aes(ymin = meanScore,ymax=meanScore, col=Sample), width = 0.5)+
#   theme_pubr() +
#   theme(aspect.ratio = 1)
# p

# There may a difference between the cohorts when taking all clonotypes, but this test is not stringent enough
rem_meanScores <- filter(neotb_grouped, survival == "Non-relapsed")$meanScore
rel_meanScores <- filter(neotb_grouped, survival == "Relapsed")$meanScore
# T-test
t_test_result <- t.test(rem_meanScores, rel_meanScores)
# Wilcox test (more commonly used in the literature)
w_test_result <- wilcox.test(rem_meanScores, rel_meanScores)
# Fold change
fc <- median(rem_meanScores)/median(rel_meanScores)

# # Sina plot
# q <- ggplot(neotb_grouped, aes(x=sample_id, y=meanScore)) +
#   geom_sina(aes(size = prop, color = sample_id, group = sample_id), scale = "width")+
#   geom_violin(alpha=0, scale = "width", draw_quantiles = 0.5) +
#   theme_pubr() +
#   theme(aspect.ratio = 1,
#         axis.text.x = element_text(angle = 45, hjust = 1))
# 
# q

# Sina plot, grouped by cohort
r <- ggplot(neotb_grouped, aes(x=survival, y=meanScore)) +
  geom_sina(aes(size = prop, color = sample_id, group = survival), scale = "width")+
  geom_violin(alpha=0, scale = "width", draw_quantiles = 0.5) +
  #ggtitle(label = unique(neoantigen$celltype)) + # this should be changed in the final script
  theme_pubr() +
  theme(aspect.ratio = 1)+
  annotate("text", x = 1.5, y = max(neotb_grouped$meanScore)-0.02,
           label = paste0("Fold change: ", round(fc, 4), "\n",
             "p = ", signif(w_test_result$p.value, digits = 4)), size=4.5)

r

pdf("Neocd8cells.pdf", width = 10, height = 8)
r
dev.off()




