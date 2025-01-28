# Nurefsan Sariipek
# Date: January 22nd, 2024
# Updated at january 28th,2025
# Analyzing post 3-6 months samples and subsetting donor/host CD4/CD8 compartments using subsampling based on cell numbers which is different than scRepertoire built in subsampling which can be found on 6.1 script
# Updated at 240418 by NS

# Load the libraries
library(scRepertoire)
library(Seurat)
library(ggpubr)
library(tidyverse)
library(janitor)
library(rstatix)
library(Hmisc)
library(glue)
#library(googleCloudStorageR)
library(RColorBrewer)

# Empty environment
rm(list=ls())

# # Set working directory 
# setwd("~/TP53_ImmuneEscape/6_TCR_Diversity/")

# Parameters to interact with Google bucket, this part only needed for Terra
gcs_global_bucket("fc-3783b423-62ac-4c69-8c2f-98cb0ee4503b")
# Check if you can list the objects. In Terra, you may need to authenticate using gcs_auth(). In VM, this did not work - hence the alternative function on line 65.
gcs_list_objects()

# Define samples
Samples <- c("P1665_MIX", "P1671_MIX", "P1745_MNC", "P1762_MIX", "P1764_MIX", 
             "P1804_MNC", "P1817_MIX", "P1964_MNC", "P2220_MNC", "P2332_MNC", "P2408_MNC", 
             "P2434_MNC", "P2448_MNC", "P2517_MIX", "P2518_CD3", "P2518_MNC", "P2599_CD3", "P2599_MNC", 
             "P2698_MIX", "P2745_MNC", "P2791_MNC", "P2820_MIX", "P2961_MNC", "P2977_MIX", "P2986_MNC", 
             "P2988_MNC", "P3000_MIX", "P6174_CD3", "P6174_MNC", "P25802_CD3", "P25802_MNC", "P25809_MNC"
)

  # Temporary directory to save downloaded files
  tmp_dir <- "/home/rstudio/tmp"
  dir.create(tmp_dir)
  
  # Track successfully loaded samples
  successful_samples <- c()
  
  # Process each sample and assign to dynamically named variables
  for (Sample in Samples) {
    print(Sample) # Log the sample being processed
    
    # Define file paths
    sample_path <- paste0(Sample, "/vdj-t/")
    file_name <- paste0(Sample, "_filtered_contig_annotations.csv")
    save_path <- file.path(tmp_dir, file_name)
    
    # Download the file and create variables
    tryCatch({
      gcs_get_object(
        object_name = paste0(sample_path, "filtered_contig_annotations.csv"),
        saveToDisk = save_path
      )
      
      # Read the downloaded file into a variable named after the sample
      assign(Sample, read.csv(save_path), envir = globalenv())
      successful_samples <- c(successful_samples, Sample) # Track success
    }, error = function(e) {
      message(glue("Error processing sample: {Sample}"))
    })
  }
  
  # Dynamically construct the contig_list for successfully loaded samples
  contig_list <- mget(successful_samples, envir = globalenv())
  
  # Combine all the samples into a single object
  combined <- combineTCR(
    contig_list,
    samples = successful_samples
  )
  
# Add variables. For scRepertoire below v2.0, replace variable.name with name
combined <- addVariable(combined, variable.name = "patient_id",
                        variables = c("P13", "P23", "P14", "P18", "P28", "P29", "P15", "P30", "P31", "P16", "P06", "P32", "P24", "P07", "P07", "P04", "P04", "P19", "P33", "P20", "P25", "P26", "P21", "P22", "P17", "P27", "P08", "P08", "P01", "P01", "P05"))

# Optional: merge data if the same sample was analyzed as both MNC and sorted T cells
combined2 <- do.call(rbind, combined)
combined <- split(combined2, f = combined2$ptnumber)

# Load the Seurat object and subset for T cells
# seu_merge <- readRDS("~/250128_seurat_annotated_final.rds")
# 
# # Keep only annotated T cell clusters (remove NK cells)
# t_cell_types <- c("CD4 Memory", "CD8 Memory", "CD4 Naïve", "Treg", "CD8 Effector","CD4 Effector Memory", "γδ T", "CD8 Exhausted", "CD8 Naïve")
# Tcells <- subset(x = seu_merge, subset = celltype %in% t_cell_types) 
# 
# # Add a new column for barcodes (rownames of the metadata)
# Tcells@meta.data$barcode <- rownames(Tcells@meta.data)
# 
# # Reset rownames to avoid duplication
# rownames(Tcells@meta.data) <- NULL
# 
# saveRDS(Tcells, "~/250128_Tcell_subset.rds")

#Next time just load the T cells
Tcells <- readRDS("~/250128_Tcell_subset.rds")

# Subset to keep only 3-6M remission samples from seurat object(might be unneccessary)
Tcells <- subset(x = Tcells, subset = timepoint %in% c("3","5","6"))

# Turn to a dataframe and keep only needed variables
meta = Tcells@meta.data
meta = meta %>% select(barcode, celltype, cohort, sample_status, orig.ident, sample_id, patient_id,timepoint, library_type)
rownames(meta) <- NULL
meta = meta %>% drop_na()

########### Optional subsetting for exploring different cell types #############

## Only subset CD8+ cells
# meta <- subset(x = meta, subset = celltype %in% c("CD8 Effector","CD8 Memory","CD8 Naïve","CD8 Terminally Exhausted"))
# meta <- subset(x = meta, subset = celltype %in% c("CD8 Effector"))
## Only subset CD4+ cells
# meta <- subset(x = meta, subset = celltype %in% c("CD4 Memory","CD4 Naïve","Treg"))

####Add Souporcell information####
# # Load the metadata that contains souporcell information
# combined_df <- read_csv(paste0(my_wd, "TP53_ImmuneEscape/5_Souporcell/outputs/cohort1-2_souporcell.csv"))
# # Wrangle the metadata
# combined_df$cell = gsub("_.*","", combined_df$cell)
# combined_df$id = gsub("\\.","_",combined_df$id )
# combined_df$barcode= paste0(combined_df$id, "_",combined_df$cell)
# # For this code only just select the remission samples
# combined_df = subset(x=combined_df, subset = status == "remission")
# combined_df = combined_df %>% select(barcode, assignment, orig.ident)
# # Join 2 dataframe
# meta <- left_join(meta, combined_df, by="barcode") %>% drop_na()
# 
# ##Subset host cells only####
# meta = subset(x=meta, subset = assignment == "host")
# # Remove P07 since it has 0 host cells
# meta = subset(x=meta, subset = patient_identity %in% c("pt01","pt02","pt04","pt05","pt06","pt08"))
# # Remove also from combined
# combined <- combined[c("P01_1Rem","P02_1Rem", "P04_1Rem", "P05_1Rem","P06_1Rem", "P08_1Rem")]

####Subset donor cells only####
 # meta = subset(x=meta, subset = assignment == "donor")

################  End of selection ###################

# Combine VDJ libraries "combined", and metadata from single cell object "meta" 
combined.sc = list()

for (i in names(combined)) {
  combined.sc[[i]] = combined[[i]] %>% 
    left_join(meta, by="barcode") %>% 
    drop_na(celltype)
}


View(combined.sc)
# Remove samples P20, P26 and P32 from the list since they had less than <400 TCR calls
samples_to_remove <- c("P20","P26", "P32")
combined.sc <- combined.sc[!names(combined.sc) %in% samples_to_remove]

# Verify the remaining samples
print(names(combined.sc))

# For some reason if you need the export the combined as a seurat object use the lines down below
# #turn combined.sc to a seurat object
# x <- do.call(rbind, combined.sc)
# row.names(x) <- x$barcode
# Tcells_TCR <- AddMetaData(Tcells, select(x, CTstrict))
# Save this combinedseurat object to use in other purposes.
#saveRDS(Tcells_TCR, file = "~/Tcells_tcr.rds"))

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

# Define the function that calculate each of these scores in a list (credits: Ksenia Safina)
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

# First, calculate and visualize how many calls we will exclude by subsetting the same number of cells from each sample
print(paste0("Note you are excluding ",
             round((1 - min(sapply(combined.sc, nrow)) * length(combined.sc) / sum(sapply(combined.sc, nrow))) * 100, 2),
             "% of the cells by subsetting to ", min(sapply(combined.sc, nrow)),
             " cells for each of the ", length(combined.sc), " samples."))
plot(sapply(combined.sc, nrow), pch = 16, ylim = c(0, max(sapply(combined.sc, nrow))))

# Visualize 
abline(h = min(sapply(combined.sc, nrow)), col = "red")

# Calculate the diversity with the function, keep in mind this only works with lists 
m = compute_diversity(combined.sc, "CTstrict", 1000)
View(m)
# For power calculation:
mean(m$inv.simpson[1:4]); sd(m$inv.simpson[1:4])
mean(m$inv.simpson[5:8]); sd(m$inv.simpson[5:8])


# Add more information to table to make more annotated plots
joined_tibble <- as_tibble(meta) %>% 
  select(sample_id, cohort, timepoint,sample_status , patient_id) %>% unique() %>%
  right_join(m, by = "patient_id")

# Determine y axis for visualization purposes
y_lim <- c(0,max(joined_tibble$inv.simpson))

# create a color patern, this was for discovery cohort

# mycolors <- c("maroon2","maroon2","royalblue1","darkgreen",
#               "lightgoldenrod3","darkturquoise",
#               "lightgreen","red3")

# Visualize the diversities
# initial box plots
p1 <- joined_tibble %>% filter(!duplicated(patient_id)) %>%
  mutate(cohort = gsub("cohort1", "Non-relapsed", gsub("cohort2", "Relapsed", cohort))) %>%
  ggplot(aes(x = cohort, y = inv.simpson)) +
  geom_boxplot()+
  geom_jitter(aes(color=patient_id), size = 4) +
  #scale_color_manual(values = mycolors)+
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
        stat_compare_means(aes(group = cohort), method = "wilcox.test", method.args = list(var.equal = T),
                           label = "p.format", label.x = 1.5, label.y=150, tip.length = 1, size = 6)

p1
# new bar graphs 240418
p2 <- joined_tibble %>% filter(!duplicated(ptnumber)) %>%
  mutate(cohort = gsub("cohort1", "Non-relapsed", gsub("cohort2", "Relapsed", cohort))) %>%
  ggplot(aes(x = cohort, y = inv.simpson)) + 
  geom_bar(stat="summary", fun=mean, aes(fill= cohort))+
  scale_fill_manual(values=c("skyblue1", "salmon"))+
  geom_jitter(aes(color=ptnumber), size = 4) +
  scale_color_manual(values = mycolors)+
  stat_summary(fun.data=mean_se, geom="errorbar", width=.5, linewidth=1) +
  theme_pubr() +
  coord_cartesian(ylim = y_lim) +
  ylab("Inverse Simpson Index") +
  theme(strip.text = element_text(size = 20 , color = "black", face="bold")) +
  theme(aspect.ratio = 2, axis.text.x = element_text(angle = 45, vjust= 1, hjust = 1,    size = 24, color = "black"),
        axis.title.x = element_blank(), 
        axis.text.y = element_text(size = 24),
        axis.title.y = element_text(size = 24, color = "black"),
        legend.key.size = unit(10,"mm"),
        legend.position = "none",
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20)) +
  stat_compare_means(aes(group = cohort), method = "t.test", #method.args = list(var.equal = T),
                     label = "p.format", label.x = 1.7, label.y=100, tip.length = 1, size = 6)

# Check the plot
p1
p2

# Save as a pdf
pdf("Post-transplant_all_wilcox.pdf", width = 6, height = 8)
p1
dev.off()




