
#January 17th, 2024
#Nurefsan will be using Ksenia's code for subsampling in this script everything else identical to 6.3_Allsamples.R

#Load the libraries
library(scRepertoire)
library(Seurat)
library(randomcoloR)
library(RColorBrewer)
library(ggpubr)
library(tidyverse)
library(janitor)

# Empty environment
rm(list=ls())

# For Nurefsan:
my_wd <- "/Users/dz855/Dropbox (Partners HealthCare)/ImmuneEscapeTP53/"

### Load the VDJ libraries ###

# Cohort 1
P01_0pre <- read.csv(paste0(my_wd, "Single Cell Data/2446_MNC/vdj_t/filtered_contig_annotations.csv"))

P01_1Rem <- read.csv(paste0(my_wd, "Single Cell Data/25802_MNC/vdj_t/filtered_contig_annotations.csv"))
P01_1RemT <- read.csv(paste0(my_wd, "Single Cell Data/25802_CD3/vdj_t/filtered_contig_annotations.csv"))
P01_2Rem <- read.csv(paste0(my_wd, "Single Cell Data/2645_MNC/vdj_t/filtered_contig_annotations.csv"))

P02_0pre <-  read.csv(paste0(my_wd, "Single Cell Data/1972_MNC/vdj_t/filtered_contig_annotations.csv"))
P02_0preT <- read.csv(paste0(my_wd, "Single Cell Data/1972_CD3/vdj_t/filtered_contig_annotations.csv"))
P02_1Rem <-  read.csv(paste0(my_wd, "Single Cell Data/2220_MNC/vdj_t/filtered_contig_annotations.csv"))
P02_2Rem <-  read.csv(paste0(my_wd, "Single Cell Data/2621_MNC/vdj_t/filtered_contig_annotations.csv"))
P02_2RemT <- read.csv(paste0(my_wd, "Single Cell Data/2621_CD3/vdj_t/filtered_contig_annotations.csv"))

P03_1Rem <-    read.csv(paste0(my_wd, "Single Cell Data/9185_MNC/vdj_t/filtered_contig_annotations.csv"))
P03_1RemT <-  read.csv(paste0(my_wd, "Single Cell Data/9185_CD3/vdj_t/filtered_contig_annotations.csv"))

P04_1Rem <-    read.csv(paste0(my_wd, "Single Cell Data/2599_MNC/vdj_t/filtered_contig_annotations.csv"))
P04_1RemT <-  read.csv(paste0(my_wd, "Single Cell Data/2599_CD3/vdj_t/filtered_contig_annotations.csv"))

# Cohort 2

P05_0pre <-  read.csv(paste0(my_wd, "Single Cell Data/9596_MNC/vdj_t/filtered_contig_annotations.csv"))
P05_0preT <-  read.csv(paste0(my_wd, "Single Cell Data/9596_CD3/vdj_t/filtered_contig_annotations.csv"))

P05_1Rem <-  read.csv(paste0(my_wd, "Single Cell Data/25809_MNC/vdj_t/filtered_contig_annotations.csv"))

P05_Rel <-  read.csv(paste0(my_wd, "Single Cell Data/2737_MNC/vdj_t/filtered_contig_annotations.csv"))
P05_RelT <-  read.csv(paste0(my_wd, "Single Cell Data/2737_CD3/vdj_t/filtered_contig_annotations.csv"))

P06_0pre <-   read.csv(paste0(my_wd, "Single Cell Data/2379_MNC/vdj_t/filtered_contig_annotations.csv"))
P06_0preT <- read.csv(paste0(my_wd, "Single Cell Data/2379_CD3/vdj_t/filtered_contig_annotations.csv"))
P06_1Rem <-   read.csv(paste0(my_wd, "Single Cell Data/2434_MNC/vdj_t/filtered_contig_annotations.csv"))

P07_1Rem <-  read.csv(paste0(my_wd, "Single Cell Data/2518_MNC/vdj_t/filtered_contig_annotations.csv"))
P07_1RemT <-  read.csv(paste0(my_wd, "Single Cell Data/2518_CD3/vdj_t/filtered_contig_annotations.csv"))
#due to the low cell number we excluded this sample
#P08_0pre <-  read.csv(paste0(my_wd, "Single Cell Data/4618_MNC/vdj_t/filtered_contig_annotations.csv"))

P08_1Rem <-  read.csv(paste0(my_wd, "Single Cell Data/6174_MNC/vdj_t/filtered_contig_annotations.csv"))
P08_1RemT <- read.csv(paste0(my_wd, "Single Cell Data/6174_CD3/vdj_t/filtered_contig_annotations.csv"))

P08_2Rem <-  read.csv(paste0(my_wd, "Single Cell Data/9931_MNC/vdj_t/filtered_contig_annotations.csv"))
P08_2RemT <-  read.csv(paste0(my_wd, "Single Cell Data/9931_CD3/vdj_t/filtered_contig_annotations.csv"))

P08_0Rel <-  read.csv(paste0(my_wd, "Single Cell Data/1953_MNC/vdj_t/filtered_contig_annotations.csv"))
P08_0RelT <-  read.csv(paste0(my_wd, "Single Cell Data/1953_CD3/vdj_t/filtered_contig_annotations.csv"))

# Make a list 
contig_list <- list(P01_0pre, P01_1Rem, P01_1RemT, P01_2Rem, P02_0pre, P02_0preT, P02_1Rem, P02_2Rem, P02_2RemT, P03_1Rem, P03_1RemT, P04_1Rem, P04_1RemT, P05_0pre, P05_0preT, P05_1Rem, P05_Rel, P05_RelT, P06_0pre, P06_0preT, P06_1Rem, P07_1Rem, P07_1RemT, P08_0Rel, P08_0RelT, P08_1Rem, P08_1RemT, P08_2Rem, P08_2RemT)

#combine the libraries
combined <- combineTCR(contig_list, 
                       samples = c("P01_0pre", "P01_1Rem", "P01_1RemT", "P01_2Rem", "P02_0pre", "P02_0preT", "P02_1Rem", "P02_2Rem", "P02_2RemT", "P03_1Rem", "P03_1RemT", "P04_1Rem","P04_1RemT", "P05_0pre", "P05_0preT", "P05_1Rem", "P05_Rel", "P05_RelT", "P06_0pre", "P06_0preT", "P06_1Rem", "P07_1Rem", "P07_1RemT", "P08_0Rel", "P08_0RelT", "P08_1Rem", "P08_1RemT", "P08_2Rem", "P08_2RemT"))   

#add this variable to combine MNC and T selected libraries later in the script
combined <- addVariable(combined, variable.name = "ptnumber",
                        variables = c("P01-0","P01-1","P01-1","P01-2","P02-0","P02-0","P02-1","P02-2","P02-2","P03-1","P03-1","P04-1","P04-1","P05-0","P05-0","P05-1","P05-2","P05-2", "P06-0","P06-0", "P06-1","P07-1","P07-1","P08-3","P08-3","P08-1","P08-1","P08-2","P08-2"))

combined <- addVariable(combined, variable.name = "cohort",
                        variables = c("cohort1","cohort1","cohort1","cohort1","cohort1","cohort1","cohort1","cohort1","cohort1","cohort1","cohort1","cohort1","cohort1","cohort2","cohort2","cohort2","cohort2","cohort2","cohort2","cohort2","cohort2","cohort2","cohort2","cohort2","cohort2","cohort2","cohort2","cohort2","cohort2")) 

# Load the Seurat object subsetted for T cells
Tcells <- readRDS(paste0(my_wd, "AnalysisNurefsan/RDS files/Tcellsfinal.rds"))

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
meta = meta %>% mutate(barcode = paste0(fullbc)) %>%
dplyr::select(barcode, celltype, cohort, orig.ident, id, Sample, groups, patient_identity, timepoint)

  
# Combine VDJ libraries "combined", and metadata from single cell object "meta" 
combined.sc = list()
  
for (i in names(combined)) {
    combined.sc[[i]] = combined[[i]] %>% 
    left_join(meta, by="barcode") %>% 
    drop_na(celltype)
    }

View(combined.sc)

#this section down below is from sc-repertoire inside of clonal diversity function, different diversity indices
.diversityCall <- function(data) {
  shannon <- .shannon(data[,"Freq"])
  inv_simpson <- .invsimpson(data[,"Freq"])
  norm_entropy <- .normentropy(data[,"Freq"]) 
  gini_simpson <- .ginisimpson(data[,"Freq"]) 
  chao1 <- .chao1(data[,"Freq"])
  ACE <- .ACE(data[,"Freq"])
  out <- c(shannon, inv_simpson, norm_entropy, gini_simpson, chao1,ACE)
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

# Implement the function that Ksenia computed
compute_diversity = function(df_list, cloneCall, n.boots) {
  mat = NULL
  min_n = min(sapply(df_list, nrow))
  for (i in seq_along(df_list)) {
    data = df_list[[i]]
    mat_a = NULL
    for (j in seq(seq_len(n.boots))) {
      x = sample_n(data, min_n)
      x = as.data.frame(table(x[,cloneCall]))
      sample = .diversityCall(x)
      mat_a = rbind(mat_a, sample)
    }
    mat_a[is.na(mat_a)] = 0
    mat_b = colMeans(mat_a)
    mat_b = as.data.frame(t(mat_b))
    mat = rbind(mat, mat_b)
  }
  colnames(mat) = c("shannon", "inv.simpson", "norm.entropy", "gini.simpson", "chao1", "ACE")
  mat$id = names(df_list)
                    
  return(mat)
}

#calculate the diversity with the function, keep in mind only works with lists
m = compute_diversity(combined.sc,"CTstrict", 1000)
View(m)

#add more information to calculation to make more annotated plots
m = m %>% left_join(meta, by="id") 

#Visualization of inverse simpson index longitudinally 
order =c("pre-transplant", "postTx_3-6m", "rem>6m", "relapse")  
m$groups= factor(m$groups, levels = order)

#timepoint is actually numerical number of the sample's timeline post-transplant another option to visualize the samples longitudinally but we have to change it to numerical using the code down below
m$timepoint = as.numeric(levels(m$timepoint)) [m$timepoint]

m1 = subset(x=m, subset = cohort == "cohort1")
m2 =subset(x=m, subset = cohort == "cohort2")



# Cohort 1
p1 <- ggplot(data = m1, aes(x = timepoint, y = inv.simpson, group = "id")) + geom_line(aes(color =id)) +facet_grid(. ~ patient_identity) + geom_point() +
  theme(aspect.ratio = 1, axis.text.x = element_text(angle = 45, vjust= 1, hjust = 1, size = 10, color = "black"),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 15, color = "black"),
        legend.key.size = unit(2,"mm"),
        legend.position = "right",
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8))    
p1
# Cohort 2
p2 <- ggplot(data = m2, aes(x = groups, y = inv.simpson, group = "id")) + geom_line(aes(color =id)) +facet_grid(. ~ patient_identity) + geom_point() +
  theme(aspect.ratio = 1, axis.text.x = element_text(angle = 45, vjust= 1, hjust = 1, size = 10, color = "black"),
      axis.title.x = element_blank(),
      axis.text.y = element_text(size = 15),
      axis.title.y = element_text(size = 15, color = "black"),
      legend.key.size = unit(2,"mm"),
      legend.position = "right",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 8)) +
theme_bw()

p2



