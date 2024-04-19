# Nurefsan Sariipek, 
# Date January 22nd, 2024
# Analyzing post 3-6 months samples and subsetting donor/host CD4/CD8 compartments using subsampling based on cell numbers which is different than scRepertoire built in subsampling which can be found on 6.1 script

# Load the libraries
library(scRepertoire)
library(Seurat)
library(ggpubr)
library(tidyverse)
library(janitor)
library(rstatix)

# Empty environment
rm(list=ls())

# For Nurefsan:
my_wd <- "/Users/dz855/Dropbox (Partners HealthCare)/ImmuneEscapeTP53/"
# For Peter:
#my_wd <- "~/DropboxMGB/Projects/ImmuneEscapeTP53/"
# Load the csv files from cohort 1 and 2 from 3-6mo remission period + from cohort 3 remission period (2 patient only)

P01_1Rem <- read.csv(paste0(my_wd, "Single Cell Data/25802_MNC/vdj_t/filtered_contig_annotations.csv"))
P01_1RemT <- read.csv(paste0(my_wd, "Single Cell Data/25802_CD3/vdj_t/filtered_contig_annotations.csv"))
P01_2Rem <- read.csv(paste0(my_wd, "Single Cell Data/2645_MNC/vdj_t/filtered_contig_annotations.csv"))
P02_1Rem <- read.csv(paste0(my_wd, "Single Cell Data/2220_MNC/vdj_t/filtered_contig_annotations.csv"))
P04_1Rem <- read.csv(paste0(my_wd, "Single Cell Data/2599_MNC/vdj_t/filtered_contig_annotations.csv"))
P04_1RemT <- read.csv(paste0(my_wd, "Single Cell Data/2599_CD3/vdj_t/filtered_contig_annotations.csv"))
P05_1Rem <- read.csv(paste0(my_wd, "Single Cell Data/25809_MNC/vdj_t/filtered_contig_annotations.csv"))
P06_1Rem <- read.csv(paste0(my_wd, "Single Cell Data/2434_MNC/vdj_t/filtered_contig_annotations.csv"))
P07_1Rem <- read.csv(paste0(my_wd, "Single Cell Data/2518_MNC/vdj_t/filtered_contig_annotations.csv"))
P07_1RemT <- read.csv(paste0(my_wd, "Single Cell Data/2518_CD3/vdj_t/filtered_contig_annotations.csv"))
P08_1Rem <- read.csv(paste0(my_wd, "Single Cell Data/6174_MNC/vdj_t/filtered_contig_annotations.csv"))
P08_1RemT <- read.csv(paste0(my_wd, "Single Cell Data/6174_CD3/vdj_t/filtered_contig_annotations.csv"))
P09_1Rem <-  read.csv(paste0(my_wd, "Single Cell Data/1732_MNC/vdj_t/filtered_contig_annotations.csv"))
P11_1Rem <- read.csv(paste0(my_wd, "Single Cell Data/6244_MNC/vdj_t/filtered_contig_annotations.csv"))
P11_1RemT <- read.csv(paste0(my_wd, "Single Cell Data/6244_CD3/vdj_t/filtered_contig_annotations.csv"))

# Make a list 
contig_list <- list(P01_1Rem, P01_1RemT, 
                    P01_2Rem, P02_1Rem, 
                    P04_1Rem, P04_1RemT,P05_1Rem, 
                    P06_1Rem, 
                    P07_1Rem, P07_1RemT, P08_1Rem, P08_1RemT,
                    P09_1Rem,P11_1Rem, P11_1RemT)

# Combine all the samples
combined <- combineTCR(contig_list,
                       samples = c("P01_1Rem", "P01_1RemT", "P01_2Rem", "P02_1Rem", 
                                   "P04_1Rem", "P04_1RemT","P05_1Rem",
                                   "P06_1Rem", 
                                   "P07_1Rem", "P07_1RemT", "P08_1Rem", "P08_1RemT", "P09_1Rem",  "P11_1Rem","P11_1RemT"))


# Add variables
combined <- addVariable(combined, name = "ptnumber",
                        variables = c("P01_1Rem", "P01_1Rem",
                                      "P01_2Rem", "P02_1Rem", 
                                      "P04_1Rem", "P04_1Rem", "P05_1Rem", "P06_1Rem",  "P07_1Rem", "P07_1Rem", "P08_1Rem", "P08_1Rem", "P09_1Rem", "P11_1Rem", "P11_1Rem"))

# Add variables
combined <- addVariable(combined, name = "cohort",
                        variables = c("cohort1","cohort1","cohort1","cohort1","cohort1","cohort1",
                                      "cohort2","cohort2","cohort2","cohort2","cohort2","cohort2", "cohort3", "cohort3", "cohort3"))


# Optional: merge data if the same sample was analyzed as both MNC and sorted T cells
combined2 <- do.call(rbind, combined)
combined <- split(combined2, f = combined2$ptnumber)

# Load the Seurat object subsetted for T cells
Tcells <- readRDS(paste0(my_wd, "AnalysisNurefsan/RDS files/Tcellsfinal.rds"))

## UMAP dimensions are lost in the previous one check this one
#Tcells <- readRDS(paste0(my_wd, "Trash/old RDS files/Tcellsubset.rds"))

# Keep only annotated T cell clusters (remove NK cells)
Tcells <- subset(x = Tcells, subset = seurat_clusters %in% c(0,1,2,3,4,5,6,7,9,10,11,12,14)) 

# Wrangle the metadata to make it compatible with the TCR metadata (combined)
Tcells_meta <- Tcells@meta.data
Tcells_meta$renamedcells <- gsub("-\\d+_\\d+","", rownames(Tcells@meta.data))
Tcells_meta <- Tcells_meta %>% mutate(fullbc = paste0(Tcells_meta$id, "_", renamedcells, "-1"))
Tcells@meta.data <- Tcells_meta

# Find duplicated barcodes and drop from Seurat object
dup_cells <- Tcells_meta$fullbc[duplicated(Tcells_meta$fullbc)]
UniqueBCs <- Tcells_meta[! Tcells_meta$fullbc %in% dup_cells,]$fullbc 
Tcells <- subset(x = Tcells, subset = fullbc %in% UniqueBCs)
Tcells <- RenameCells(Tcells, new.names = UniqueBCs)

# Turn to a dataframe and keep only needed variables
meta = Tcells@meta.data
meta = meta %>% mutate(barcode = paste0(fullbc), ptnumber = paste0(Sample))
meta = meta %>% select(barcode, celltype, cohort, orig.ident, id, ptnumber, groups, patient_identity, timepoint)
rownames(meta) <- NULL
meta = meta %>% drop_na()

# Optional subsetting for exploring different cell types #
####Subset celltypes#####
# Only subset CD8+ cells
meta <- subset(x = meta, subset = celltype %in% c("CD8 Effector","CD8 Memory","CD8 Naïve","CD8 Terminally Exhausted"))
#each CD8 celltype
meta <- subset(x = meta, subset = celltype == "CD8 Memory")
meta <- subset(x = meta, subset = celltype == "CD8 Naïve")
meta <- subset(x = meta, subset = celltype == "CD8 Effector")
meta <- subset(x = meta, subset = celltype == "CD8 Terminally Exhausted")

# Only subset CD4+ cells
meta <- subset(x = meta, subset = celltype %in% c("CD4 Memory","CD4 Naïve","Treg"))

###Add Souporcell information####
# Load the metadata that contains souporcell information
combined_df <- read_csv(paste0(my_wd, "AnalysisNurefsan/Souporcell/outputs/cohort1-2_souporcell.csv"))
# Wrangle the metadata 
combined_df$cell = gsub("_.*","", combined_df$cell)
combined_df$id = gsub("\\.","_",combined_df$id )
combined_df$barcode= paste0(combined_df$id, "_",combined_df$cell)
# For this code only just select the remission samples
combined_df = subset(x=combined_df, subset = status == "remission")
combined_df = combined_df %>% select(barcode, assignment, orig.ident)
# Join 2 dataframe
meta <- left_join(meta,combined_df, by="barcode") %>% drop_na()

###Subset host cells only####
meta = subset(x=meta, subset = assignment == "host")
# Remove P07 since it has 0 host cells
meta = subset(x=meta, subset = patient_identity %in% c("pt01","pt02","pt04","pt05","pt06","pt08"))
# Remove also from combined
combined <- subsetClones(combined, 
                         name = "ptnumber", 
                         variables = c("P01_1Rem","P02_1Rem", "P04_1Rem",  "P05_1Rem","P06_1Rem","P07_1Rem", "P08_1Rem"))

###Subset donor cells only####
meta = subset(x=meta, subset = assignment == "donor")

####End of selection####
# Combine VDJ libraries "combined", and metadata from single cell object "meta" 
combined.sc = list()

for (i in names(combined)) {
  combined.sc[[i]] = combined[[i]] %>% 
    left_join(meta, by="barcode") %>% 
    drop_na(celltype)
}

View(combined.sc)

# This section is from the scRepertoire clonal diversity function. It calculates different diversity indices.
.diversityCall <- function(data) {
  shannon <- .shannon(data[,"Freq"])
  inv_simpson <- .invsimpson(data[,"Freq"])
  norm_entropy <- .normentropy(data[,"Freq"]) 
  gini_simpson <- .ginisimpson(data[,"Freq"]) 
  chao1 <- .chao1(data[,"Freq"])
  ACE <- .ACE(data[,"Freq"])
  out <- c(shannon, inv_simpson, norm_entropy, gini_simpson, chao1, ACE)
  return(out)
}

.shannon <- function(p){
  p <- p[which(p > 0)]
  p <- p / sum(p)
  p <- p[which(p > 0)]
  return(-sum(p * log(p)))
}
.normentropy <- function(p){
  p <- p[which(p > 0)]
  p <- p / sum(p)
  p <- p[which(p > 0)]
  return(-sum(p * log(p)) / log(length(p)))
}
.invsimpson <- function(p){
  p <- p[which(p > 0)]
  p <- p / sum(p)
  p <- p[which(p > 0)]
  return(1 / sum(p^2))
}
.ginisimpson <- function(p){
  p <- p[which(p > 0)]
  p <- p / sum(p)
  p <- p[which(p > 0)]
  return(1 - sum(p^2))
}

.chao1 <- function(p){
  n1 <- sum(p == 1)
  n2 <- sum(p == 2)
  S_obs <- length(p)
  # Chao1 index calculation
  if(n1 > 1 && n2 > 0) {
    chao1 <- S_obs + (n1 * (n1 - 1)) / (2 * (n2 + 1))
  } else {
    # In cases where n1 <= 1 or n2 == 0, Chao1 is undefined
    chao1 <- NA
  }
  return(chao1)
}

.ACE <- function(p) {
  q <- 10
  S_abund <- sum(p > q)
  rare_data <- p[p <= q]
  S_rare <- length(rare_data)
  n_rare <- sum(rare_data)
  
  # Calculate C_ACE
  C_ACE <- sum(p) / n_rare
  
  # Calculate gamma
  gamma <- 0
  for(i in seq_len(q)) {
    f_i <- sum(rare_data == i)
    gamma <- gamma + (1 - i / q)^f_i
  }
  
  # Calculate ACE
  ACE <- S_abund + (S_rare / C_ACE) + (1 - C_ACE) * gamma
  return(ACE)
}

# Define the function that Ksenia wrote
compute_diversity = function(df_list, cloneCall, n.boots) {
  # df_list <- combined.sc
  # cloneCall <- "CTstrict"
  # n.boots <- 1000
  mat = NULL
  min_n = min(sapply(df_list, nrow))
  for (i in seq_along(df_list)) {
    data = df_list[[i]]
    mat_a = NULL
    for (j in seq(seq_len(n.boots))) {
      x = slice_sample(data, n = min_n)
      y = as.data.frame(table(x[,cloneCall]))
      sample = .diversityCall(y)
      mat_a = rbind(mat_a, sample)
    }
    mat_a[is.na(mat_a)] = 0
    mat_b = colMeans(mat_a)
    mat_b = as.data.frame(t(mat_b))
    mat = rbind(mat, mat_b)
  }
  colnames(mat) = c("shannon", "inv.simpson", "norm.entropy", "gini.simpson", "chao1", "ACE")
  mat$ptnumber = names(df_list)
  
  return(mat)
}

# First, calculate and visualize how many cells we will exclude by taking the minimum number
print(paste0("Note you are excluding ",
             round((1 - min(sapply(combined.sc, nrow)) * length(combined.sc) / sum(sapply(combined.sc, nrow))) * 100, 2),
             "% of the cells by subsetting to ", min(sapply(combined.sc, nrow)),
             " cells for each of the ", length(combined.sc), " samples."))
plot(sapply(combined.sc, nrow), pch = 16, ylim = c(0, max(sapply(combined.sc, nrow))))
abline(h = min(sapply(combined.sc, nrow)), col = "red")
# Calculate the diversity with the function, keep in mind only works with lists 
m = compute_diversity(combined.sc, "CTstrict", 1000)
View(m)

# Add more information to calculation to make more annotated plots
joined_tibble <- as_tibble(meta) %>% 
  select(id, cohort, timepoint, groups, patient_identity, ptnumber) %>% unique() %>%
  right_join(m, by = "ptnumber")

# Determine y axis
y_lim <- c(0,max(joined_tibble$inv.simpson))

# Visualize the diversities

p1 <- joined_tibble %>% filter(!duplicated(ptnumber)) %>%
  mutate(cohort = gsub("cohort1", "Non-relapsed", gsub("cohort3", "Early-Relapsed", gsub("cohort2", "Relapsed", cohort)))) %>%
  ggplot(aes(x = cohort, y = inv.simpson)) + 
  geom_boxplot()+
  geom_jitter(aes(color=patient_identity), size = 4) +
  theme_bw() +
  coord_cartesian(ylim = y_lim) +
  ylab("Inverse Simpson Index") +
  theme(strip.text = element_text(size = 14, color = "black", face="bold")) +
  theme(aspect.ratio = 1.5, axis.text.x = element_text(angle = 45, vjust= 1, hjust = 1,    size = 15, color = "black"),
        axis.title.x = element_blank(), 
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 15, color = "black"),
        legend.key.size = unit(3,"mm"),
        legend.position = "right",
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12)) +
  stat_compare_means(aes(group = cohort), method = "t.test", method.args = list(var.equal = T),
                     label = "p.signif", label.x = 2, label.y=20, tip.length = 1, size = 6)
# Check the plot
p1
# Save as a pdf
pdf("Post-transplantcdeffector.pdf", width = 5, height = 7.5)
p1
dev.off()




