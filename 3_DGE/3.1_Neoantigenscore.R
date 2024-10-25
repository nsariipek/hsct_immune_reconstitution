## Neantigen score+TCRs
# Nurefsan Sariipek, 24-07-01

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

# For Nurefsan:
my_wd <- "/Users/dz855/Dropbox (Partners HealthCare)/ImmuneEscapeTP53/"
# For Peter
#my_wd <- "~/DropboxMGB/Projects/ImmuneEscapeTP53/TP53_ImmuneEscape/3_DGE/3.1_Neoantigenscore/"

# Load the seurat object from 6.1 script end of line 192 which has the TCR+ scRNA combined object
combined <- readRDS(paste0(my_wd, "RDS files/Combined.rds"))
# Load antigen scores from Paper(nurefsan insert the paper) as a data table 
cd4neoA <- read.csv("TP53_ImmuneEscape/3_DGE/3.1_Neoantigenscore/signatures/cd4.csv")
cd8neoA <- read.csv("TP53_ImmuneEscape/3_DGE/3.1_Neoantigenscore/signatures/cd8.csv")
neoA <- read.csv("TP53_ImmuneEscape/3_DGE/3.1_Neoantigenscore//signatures/neoantigen.csv")

# Subset for different celltypes before adding the antigen score
#cd8cells <- subset(x = combined, subset = celltype %in% c("CD8 Effector","CD8 Memory","CD8 Naïve","CD8 Terminally Exhausted"))
 cd8ef <- subset(x = combined, subset = celltype == "CD8 Effector")
# cd8mem <- subset(x = combined, subset = celltype =="CD8 Memory")
# cd8na <- subset(x = combined, subset = celltype=="CD8 Naïve")
# cd8ex <- subset(x = combined, subset = celltype=="CD8 Terminally Exhausted")
 cd4cells <- subset(x = combined, subset = celltype %in% c("CD4 Memory","CD4 Naïve","Treg"))
# treg <- subset(x = combined, subset = celltype == "Treg")
# cd4mem <- subset(x = combined, subset = celltype == "CD4 Memory")
# cd4na <- subset(x = combined, subset = celltype == "CD4 Naïve")

# Add to the seurat object
neoantigen <- AddModuleScore(object =  cd8ef,
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

# Sanity check, visualize the violin plots
neoantigen <- SetIdent(neoantigen, value = "Sample")
VlnPlot(neoantigen, features = "neoantigen1", split.by = "cohort", sort = "increasing")

# Turn into tibble
neoantigen <- as_tibble(neoantigen@meta.data, rownames = "cell") 

neotb <- neoantigen%>%
  select(Sample, CTstrict, neoantigen1)

# For each TCR, calculate relative clonotype size (grouped by sample)
neotb_grouped <- neotb %>%
  group_by(Sample) %>%
  mutate(n_total = n()) %>%
  ungroup() %>%
  group_by(Sample, CTstrict) %>%
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

# add a column for cohort information
neotb_grouped$cohort <- case_when(grepl("P01|P02|P04",neotb_grouped$Sample) ~ "cohort1", grepl("P05|P06|P07|P08",neotb_grouped$Sample) ~ "cohort2")

# There may a difference between the cohorts when taking all clonotypes, but this test is not stringent enough
cohort1_meanScores <- filter(neotb_grouped, cohort == "cohort1")$meanScore
cohort2_meanScores <- filter(neotb_grouped, cohort == "cohort2")$meanScore
# T-test
t_test_result <- t.test(cohort1_meanScores, cohort2_meanScores)
# Wilcox test (more commonly used in the literature)
w_test_result <- wilcox.test(cohort1_meanScores, cohort2_meanScores)
# Fold change
fc <- median(cohort1_meanScores)/median(cohort2_meanScores)

# Sina plot
q <- ggplot(neotb_grouped, aes(x=Sample, y=meanScore)) +
  geom_sina(aes(size = prop, color = Sample, group = Sample), scale = "width")+
  geom_violin(alpha=0, scale = "width", draw_quantiles = 0.5) +
  theme_pubr() +
  theme(aspect.ratio = 1,
        axis.text.x = element_text(angle = 45, hjust = 1))

q

# Sina plot, grouped by cohort
r <- ggplot(neotb_grouped, aes(x=cohort, y=meanScore)) +
  geom_sina(aes(size = prop, color = Sample, group = cohort), scale = "width")+
  geom_violin(alpha=0, scale = "width", draw_quantiles = 0.5) +
  ggtitle(label = unique(neoantigen$celltype)) + # this should be changed in the final script
  theme_pubr() +
  theme(aspect.ratio = 1) +
  annotate("text", x = 1.5, y = max(neotb_grouped$meanScore)-0.02,
           label = paste0("Fold change: ", round(fc, 4), "\n",
             "p = ", signif(w_test_result$p.value, digits = 4)), size=4.5)

r

pdf("cd8eff.pdf", width = 8, height = 6)
r
dev.off()




