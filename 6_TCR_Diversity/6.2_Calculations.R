#Nov,2023
# Trying to calculate indeces by ourselves with this code 
library(scRepertoire)
library(Seurat)
library(randomcoloR)
library(RColorBrewer)
library(ggpubr)
library(tidyverse)
library(janitor)
library(vegan)

# Empty environment
rm(list=ls())

# For Nurefsan:
my_wd <- "/Users/dz855/Dropbox (Partners HealthCare)/ImmuneEscapeTP53/"
# Load the Seurat object subsetted for T cells
Tcells <- readRDS(paste0(my_wd, "AnalysisNurefsan/RDS files/Tcellsfinal.rds"))


# Keep only annotated T cell clusters (remove NK cells)
Tcells <- subset(x = Tcells, subset = seurat_clusters %in% c(0,1,2,3,4,5,6,7,9,10,11,12,14)) 
Tcells_meta <- Tcells@meta.data
Tcells_meta$renamedcells <- gsub("-\\d+_\\d+","", rownames(Tcells@meta.data))
Tcells_meta <- Tcells_meta %>% mutate(fullbc = paste0(Tcells_meta$id, "_", renamedcells, "-1"))
Tcells@meta.data <- Tcells_meta

# Find duplicated barcodes and drop from Seurat object
dup_cells <- Tcells_meta$fullbc[duplicated(Tcells_meta$fullbc)]
UniqueBCs <- Tcells_meta[! Tcells_meta$fullbc %in% dup_cells,]$fullbc 
Tcells <- subset(x = Tcells, subset = fullbc %in% UniqueBCs)
Tcells <- RenameCells(Tcells, new.names = UniqueBCs)

# Combine TCR and sc-RNAseq data
Tcells_combined <- combineExpression(combined, Tcells, 
                                     cloneCall = "strict",
                                     proportion = FALSE,
                                     filterNA = T,
                                     clone=c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=500))

Tcells_combined_tib <- as_tibble(Tcells_combined@meta.data, rownames = "cell")
Tcells_combined_tib <- Tcells_combined_tib %>%
  mutate(pt_timepoint = gsub("RemT", "Rem", id)) %>%
  mutate(pt_timepoint_ct = paste0(pt_timepoint, "_", CTstrict))

# Now calculate the Frequency / clone size
Tcells_combined_tib <- Tcells_combined_tib %>% group_by(pt_timepoint_ct) %>% add_count() %>%
  ungroup() %>% arrange(pt_timepoint_ct) %>% arrange(pt_timepoint, n)


# Load the metadata that contains souporcell information
combined_df <- read_csv(paste0(my_wd, "AnalysisNurefsan/Souporcell/outputs/cohort1-2_souporcell.csv"))


# Wrangle the metadata to 
combined_df$cell = gsub("_.*","", combined_df$cell)
combined_df$id = gsub("\\.","_",combined_df$id )
combined_df$cell= paste0(combined_df$id, "_",combined_df$cell)

# Merge 2 metadata
newdf <- Tcells_combined_tib %>% 
  left_join(combined_df, by ="cell") %>% 
  drop_na()

unique(newdf$pt_timepoint)

donorcells%>% tabyl(id.x)

# Subset donor cells only

donorcells <- subset(x = newdf, subset = assignment == "donor")
hostcells <- subset(x = newdf, subset = assignment == "host")
donorcd4 <- subset(x = donorcells, subset = celltype.y %in% c("CD4 Memory","CD4 Naïve","Treg"))
hostcd4 <- subset(x = hostcells, subset = celltype.y %in% c("CD4 Memory","CD4 Naïve","Treg"))
donorcd8 <- subset(x = donorcells, subset = celltype.y %in% c("CD8 Memory","CD8 Naïve","CD8 Effector", "CD8 Terminally Exhausted"))
hostcd8 <- subset(x = hostcells, subset = celltype.y %in% c("CD8 Memory","CD8 Naïve","CD8 Effector", "CD8 Terminally Exhausted"))

newdf %>% tabyl(pt_timepoint)

# Add assignment calls to the Seurat metadata
Tcells_combined <- AddMetaData(Tcells_combined, data.frame(select(newdf, cell, assignment), row.names = "cell"))
Tcells_combined$assignment %>% tabyl

ds = data.frame()

for (i in c(1:1000)) {
Tcells_subset_tib <- newdf %>% group_by(pt_timepoint) %>% slice_sample(n =190) %>%
  select(-n) %>% group_by(pt_timepoint_ct) %>% add_count() %>%
  ungroup() %>% arrange(pt_timepoint_ct) %>% arrange(pt_timepoint, n)


df <-Tcells_subset_tib %>% 
  group_by(pt_timepoint,CTstrict) %>% 
  summarise(n= sum(n)) %>% 
  ungroup %>% 
  spread(CTstrict,n, fill=0) 


rand <- Tcells_subset_tib%>%
  uncount(clonalFrequency) %>%
  mutate(name = sample(pt_timepoint)) %>%
  count(name, name="value")
  
richness <- function(x){
  
  # r <- sum(x > 0)
  # return(r)
  
  sum(x>0)
}

shannon <- function(x){
  
  rabund <- x[x>0]/sum(x)
  -sum(rabund * log(rabund))
  
}

simpson <- function(x){
  
  n <- sum(x)
  
  # sum(x * (x-1) / (n * (n-1)))
  1 - sum((x/n)^2)
}


rand %>%
  group_by(name) %>%
  summarize(sobs = richness(value),
            simpson = simpson(value),
            invsimpson = 1/simpson)
           # n = sum(clonalFrequency)) 



d = diversity(df[,-1], "invsimpson")


ds = rbind(ds,d)
}

colMeans(ds) 

