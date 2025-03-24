#January 17th, 2024
# Nurefsan will be using Ksenia's code for subsampling in this script everything else identical to 6.3_Allsamples.R

# Load the libraries
library(scRepertoire)
library(Seurat)
library(randomcoloR)
library(RColorBrewer)
library(ggpubr)
library(tidyverse)
library(janitor)
library(cowplot)

# Empty environment
rm(list=ls())

# For Nurefsan:
my_wd <- "/Users/dz855/Dropbox (Partners HealthCare)/ImmuneEscapeTP53/"
# For Peter:
#my_wd <- "~/DropboxMGB/Projects/ImmuneEscapeTP53/"

### Load the VDJ libraries ###

# Cohort 1
P01_0pre <- read.csv(paste0(my_wd, "TCR files/2446_MNC_filtered_contig_annotations.csv"))

P01_1Rem <- read.csv(paste0(my_wd, "TCR files/25802_MNC_filtered_contig_annotations.csv"))
P01_1RemT <- read.csv(paste0(my_wd, "TCR files/25802_CD3_filtered_contig_annotations.csv"))
P01_2Rem <- read.csv(paste0(my_wd, "TCR files/2645_MNC_filtered_contig_annotations.csv"))

P02_0pre <-  read.csv(paste0(my_wd, "TCR files/1972_MNC_filtered_contig_annotations.csv"))
P02_0preT <- read.csv(paste0(my_wd, "TCR files/1972_CD3_filtered_contig_annotations.csv"))
P02_1Rem <-  read.csv(paste0(my_wd, "TCR files/2220_MNC_filtered_contig_annotations.csv"))
P02_2Rem <-  read.csv(paste0(my_wd, "TCR files/2621_MNC_filtered_contig_annotations.csv"))
P02_2RemT <- read.csv(paste0(my_wd, "TCR files/2621_CD3_filtered_contig_annotations.csv"))

P03_1Rem <-    read.csv(paste0(my_wd, "TCR files/9185_MNC_filtered_contig_annotations.csv"))
P03_1RemT <-  read.csv(paste0(my_wd, "TCR files/9185_CD3_filtered_contig_annotations.csv"))

P04_1Rem <-    read.csv(paste0(my_wd, "TCR files/2599_MNC_filtered_contig_annotations.csv"))
P04_1RemT <-  read.csv(paste0(my_wd, "TCR files/2599_CD3_filtered_contig_annotations.csv"))

# Cohort 2

P05_0pre <-  read.csv(paste0(my_wd, "TCR files/9596_MNC_filtered_contig_annotations.csv"))
P05_0preT <-  read.csv(paste0(my_wd, "TCR files/9596_CD3_filtered_contig_annotations.csv"))

P05_1Rem <-  read.csv(paste0(my_wd, "TCR files/25809_MNC_filtered_contig_annotations.csv"))

P05_Rel <-  read.csv(paste0(my_wd, "TCR files/2737_MNC_filtered_contig_annotations.csv"))
P05_RelT <-  read.csv(paste0(my_wd, "TCR files/2737_CD3_filtered_contig_annotations.csv"))

P06_0pre <-   read.csv(paste0(my_wd, "TCR files/2379_MNC_filtered_contig_annotations.csv"))
P06_0preT <- read.csv(paste0(my_wd, "TCR files/2379_CD3_filtered_contig_annotations.csv"))
P06_1Rem <-   read.csv(paste0(my_wd, "TCR files/2434_MNC_filtered_contig_annotations.csv"))

P07_1Rem <-  read.csv(paste0(my_wd, "TCR files/2518_MNC_filtered_contig_annotations.csv"))
P07_1RemT <-  read.csv(paste0(my_wd, "TCR files/2518_CD3_filtered_contig_annotations.csv"))
#due to the low cell number we excluded this sample
#P08_0pre <-  read.csv(paste0(my_wd, "Single Cell Data/4618_MNC/vdj_t/filtered_contig_annotations.csv"))

P08_1Rem <-  read.csv(paste0(my_wd, "TCR files/6174_MNC_filtered_contig_annotations.csv"))
P08_1RemT <- read.csv(paste0(my_wd, "TCR files/6174_CD3_filtered_contig_annotations.csv"))

P08_2Rem <-  read.csv(paste0(my_wd, "TCR files/9931_MNC_filtered_contig_annotations.csv"))
P08_2RemT <-  read.csv(paste0(my_wd, "TCR files/9931_CD3_filtered_contig_annotations.csv"))

P08_0Rel <-  read.csv(paste0(my_wd, "TCR files/1953_MNC_filtered_contig_annotations.csv"))
P08_0RelT <-  read.csv(paste0(my_wd, "TCR files/1953_CD3_filtered_contig_annotations.csv"))

# Make a list 
contig_list <- list(P01_0pre, P01_1Rem, P01_1RemT, P01_2Rem, P02_0pre, P02_0preT, P02_1Rem, P02_2Rem, P02_2RemT, P03_1Rem, P03_1RemT, P04_1Rem, P04_1RemT, P05_0pre, P05_0preT, P05_1Rem, P05_Rel, P05_RelT, P06_0pre, P06_0preT, P06_1Rem, P07_1Rem, P07_1RemT, P08_0Rel, P08_0RelT, P08_1Rem, P08_1RemT, P08_2Rem, P08_2RemT)

#combine the libraries
combined <- combineTCR(contig_list, 
                       samples = c("P01_0pre", "P01_1Rem", "P01_1RemT", "P01_2Rem", "P02_0pre", "P02_0preT", "P02_1Rem", "P02_2Rem", "P02_2RemT", "P03_1Rem", "P03_1RemT", "P04_1Rem","P04_1RemT", "P05_0pre", "P05_0preT", "P05_1Rem", "P05_Rel", "P05_RelT", "P06_0pre", "P06_0preT", "P06_1Rem", "P07_1Rem", "P07_1RemT", "P08_0Rel", "P08_0RelT", "P08_1Rem", "P08_1RemT", "P08_2Rem", "P08_2RemT"))   

#add this variable to combine MNC and T selected libraries later in the script
combined <- addVariable(combined, name = "ptnumber",
                        variables = c("P01_0pre","P01_1Rem","P01_1Rem", "P01_2Rem", "P02_0pre", "P02_0pre", "P02_1Rem", "P02_2Rem", "P02_2Rem", "P03_1Rem", "P03_1Rem", "P04_1Rem", "P04_1Rem", "P05_0pre", "P05_0pre", "P05_1Rem", "P05_Rel", "P05_Rel",       "P06_0pre", "P06_0pre", "P06_1Rem",   "P07_1Rem", "P07_1Rem", "P08_0Rel", "P08_0Rel", "P08_1Rem", "P08_1Rem", "P08_2Rem", "P08_2Rem"))

# Optional: merge data if the same sample was analyzed as both MNC and sorted T cells
combined2 <- do.call(rbind, combined)
combined <- split(combined2, f = combined2$ptnumber)

# Load the Seurat object subsetted for T cells
Tcells <- readRDS(paste0(my_wd, "RDS files/Tcellsfinal.rds"))

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

# Visualization of inverse simpson index longitudinally
joined_tibble$groups = factor(joined_tibble$groups, levels = c("pre-transplant", "postTx_3-6m", "rem>6m", "relapse"))

# Change timepoint from factor to numerical value in order to visualize the samples longitudinally
joined_tibble$timepoint = as.numeric(levels(joined_tibble$timepoint))[joined_tibble$timepoint]

m1 = subset(x=joined_tibble, subset = cohort == "cohort1")
m2 = subset(x=joined_tibble, subset = cohort == "cohort2")


y_lim <- c(0,max(joined_tibble$inv.simpson))
solid.line.points <- c('pre-transplant','postTx_3-6m')
dashed.line.points <- c('postTx_3-6m','postTx_3-6m','rem>6m', 'relapse')  

# Cohort 1

mycolors <- c("maroon2","royalblue1","magenta4","darkgreen")
              
p1 <- 
  ggplot(m1,aes(groups,inv.simpson, group=patient_identity)) + 
  geom_point(size =5,aes(color=patient_identity), shape=1) +
  stat_summary(data=subset(m1, groups %in% solid.line.points), fun=mean, geom="line", aes(x=groups, y=inv.simpson,linetype='solid line', color= patient_identity)) +
  stat_summary(data=subset(m1, groups %in% dashed.line.points),fun=mean, geom="line", aes(x=groups, y=inv.simpson, linetype='dashed line',color= patient_identity))+
  guides(linetype = "none")+
  scale_color_manual(values = mycolors)+
  #scale_color_brewer(palette="Set1")+
  theme_bw() +
  coord_cartesian(ylim = y_lim) +
  theme(aspect.ratio = 0.75, axis.text.x = element_text(angle = 45, vjust= 1, hjust = 1, size = 16, color = "black"),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 16),
        axis.title.y = element_text(size = 20, color = "black"),
        legend.key.size = unit(8,"mm"),
        legend.position = "right",
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16)) 


p1

# Cohort 2
mycolors <- c("lightgoldenrod3","darkturquoise",
              "lightgreen","red3")
p2 <- 
ggplot(m2,aes(groups,inv.simpson, group=patient_identity)) + 
  #facet_grid(. ~ patient_identity) + 
  geom_point(size =5,aes(color=patient_identity), shape=1) +
  stat_summary(data=subset(m2, groups %in% solid.line.points), fun=mean, geom="line", aes(x=groups, y=inv.simpson,linetype='solid line', color= patient_identity)) +
  stat_summary(data=subset(m2, groups %in% dashed.line.points),fun=mean, geom="line", aes(x=groups, y=inv.simpson, linetype='dashed line',color= patient_identity))+
  guides(linetype = "none")+
  scale_color_manual(values = mycolors)+
  #scale_color_brewer(palette="Dark2")+
theme_bw() +
  coord_cartesian(ylim = y_lim) +
  theme(aspect.ratio =0.75, axis.text.x = element_text(angle = 45, vjust= 1, hjust = 1, size = 16, color = "black"),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 16),
        axis.title.y = element_text(size = 20, color = "black"),
        legend.key.size = unit(8,"mm"),
        legend.position = "right",
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14)) 

p2
pdf("p2.pdf", width = 6, height = 4)
p2
dev.off()

pdf("p1.pdf", width = 6, height = 4)
p1
dev.off()

p <- plot_grid(p1,p2)
save_plot("longitudinaldiversity.pdf", p, ncol = 2, base_asp = 2, base_height = 4, base_width = 6)
dev.off()



