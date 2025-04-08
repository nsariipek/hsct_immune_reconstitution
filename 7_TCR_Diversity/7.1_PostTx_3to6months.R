# Nurefsan Sariipek
# Date: January 22nd, 2024
# Updated February 14, 2025
# Analyze post 3-6 months samples and subsetting donor/host CD4/CD8 compartments using subsampling based on cell numbers which is different than scRepertoire built in subsampling which can be found on 6.1 script

# Load the libraries
library(scRepertoire)
library(Seurat)
library(ggpubr)
library(tidyverse)
library(janitor)
library(rstatix)
library(Hmisc)
library(glue)
library(googleCloudStorageR) # For Terra
library(RColorBrewer)

# Empty environment
rm(list=ls())

# Set working directory 
setwd("~/TP53_ImmuneEscape/7_TCR_Diversity/")

# Setup for Peter: local. Then run line 89-117 and continue at 178
#setwd("~/DropboxMGB/Projects/ImmuneEscapeTP53/TP53_ImmuneEscape/6_TCR_Diversity")
#combined.sc <- readRDS("AuxiliaryFiles/combined.RDS")
#seu_merge <- readRDS("../AuxiliaryFiles/250128_seurat_annotated_final.rds")

# Parameters to interact with Google bucket, this part only needed for Terra
gcs_global_bucket("fc-3783b423-62ac-4c69-8c2f-98cb0ee4503b")
# Check if you can list the objects. In Terra, you may need to authenticate using gcs_auth(). In VM, this did not work - hence the alternative function on line 65.
gcs_list_objects()

# Define samples
Samples <- c("P1665_MIX", "P1671_MIX", "P1745_MNC", "P1762_MIX", "P1764_MIX", 
             "P1804_MNC", "P1817_MIX", "P1964_MNC", "P2220_MNC", "P2332_MNC", "P2408_MNC", 
             "P2434_MNC", "P2448_MNC", "P2517_MIX", "P2518_CD3", "P2518_MNC", "P2599_CD3", "P2599_MNC", "P2698_MIX", "P2745_MNC", "P2791_MNC", "P2820_MIX", "P2961_MNC", "P2977_MIX", "P2986_MNC", "P2988_MNC", "P3000_MIX", "P6174_CD3", "P6174_MNC", "P25802_CD3", "P25802_MNC","P25809_MNC")

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
                        variables = c("P13", "P23", "P14", "P18", "P28", "P29", "P15", "P30","P02", "P31", "P16", "P06", "P32", "P24", "P07", "P07", "P04", "P04", "P19", "P33", "P20", "P25", "P26", "P21", "P22", "P17", "P27", "P08", "P08", "P01", "P01", "P05"))

# Optional: merge data if the same sample was analyzed as both MNC and sorted T cells
combined2 <- do.call(rbind, combined)
combined <- split(combined2, f = combined2$patient_id)

# load the T cells
Tcells <- readRDS("~/250128_Tcell_subset.rds")

# Define TP53 MT and WT patient groups
mt_patients <- paste0("P", sprintf("%02d", c(1:12, 14, 17)))
wt_patients <- paste0("P", c(13, 15, 16, 18:33))  

Tcells@meta.data <- Tcells@meta.data %>%
  mutate(patient_id = as.factor(patient_id))%>%
  mutate(TP53 = case_when(
    patient_id %in% mt_patients ~ "MT",
    patient_id %in% wt_patients ~ "WT"))

# Subset to keep only 3-6M remission samples from seurat object(might be unneccessary)
Tcells <- subset(x = Tcells, subset = timepoint %in% c("3","5","6") & sample_status == "remission")

# Turn to a dataframe and keep only needed variables
meta = Tcells@meta.data %>% 
       select(barcode, celltype, cohort, sample_status, orig.ident, sample_id, patient_id,timepoint,survival, library_type,TP53)
rownames(meta) <- NULL
meta = meta %>% drop_na()

########### Optional subsetting for exploring different cell types #############
# # Only MT samples
# meta <- subset(x=meta, subset = TP53=="MT")
# 
# # Only subset CD8+ cells
#  meta <- subset(x = meta, subset = celltype %in% c("CD8 Memory", "CD8 Effector", "CD8 Exhausted","γδ T","CD8 Naïve"))
#  
# # Only subset CD4+ cells
# meta <- subset(x = meta, subset = celltype %in% c("Treg","CD4 Effector Memory", "CD4 Naïve","CD4 Memory"))

### Add Souporcell information ####

# # Load the metadata that contains souporcell information
# souporcell_df <- read_csv("~/final_dataset.csv")
# 
# # Wrangle the df 
# souporcell_df <- souporcell_df %>%
#   mutate(modified_barcode = paste(orig.ident, barcode, sep = "_")) %>% 
#   select(modified_barcode, origin)
# 
# # Join 2 dataframe
# meta_merged <- souporcell_df %>%
#   inner_join(meta, by = c("modified_barcode" = "barcode")) %>% drop_na()

# # Select only recipient cells
# meta <- meta_merged %>% rename(barcode=modified_barcode) %>% filter(origin == "donor")

################  End of selection ###################

# Combine VDJ libraries "combined", and metadata from single cell object "meta" 
combined.sc = list()

for (i in names(combined)) {
  combined.sc[[i]] = combined[[i]] %>% 
    left_join(meta, by = "barcode") %>% 
    drop_na(celltype) 
}

View(combined.sc)

# P06,P16,P20,P24,P26,P32
# Remove samples P20, P26 and P32 from the list since they had less than <500 TCR calls

samples_to_remove <- c("P16","P20","P26","P32","P33")
#mt_samples <- c("P01","P02","P03","P04","P05","P06","P07","P08","P10","P11","P12","P14","P17")
combined.sc <- combined.sc[!names(combined.sc) %in% samples_to_remove]
#combined.sc <- combined.sc[names(combined.sc) %in% mt_samples]

# Verify the remaining samples
print(names(combined.sc))

# You need to save RDS with TCR info for calculating neoantigen signature 
# Turn combined.sc to a seurat object
# x <- do.call(rbind, combined.sc)
# row.names(x) <- x$barcode
# Tcells_TCR <- AddMetaData(Tcells, select(x, CTstrict))
# # Save this combinedseurat object to use in other purposes.
# saveRDS(Tcells_TCR, file = "~/Tcells_TCR.rds")

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
  mat$patient_id = names(df_list)
  
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

# Add more information to table to make more annotated plots
joined_tibble <- as_tibble(meta) %>% 
  select(sample_id, cohort, timepoint,sample_status , patient_id, survival,TP53,celltype) %>% unique() %>%
  right_join(m, by = "patient_id")

# Determine y axis for visualization purposes
y_lim <- c(0,max(joined_tibble$inv.simpson))

# Visualize the diversities
# Survival colors
survival_colors <- c("Non-relapsed" = "#4775FFFF","Relapsed" = "#E64B35FF")
p2 <- joined_tibble %>%
  filter(!duplicated(patient_id)) %>%
  filter(TP53=="WT") %>%
  ggplot(aes(x = survival, y = inv.simpson, fill = survival)) +
  geom_bar(stat = "summary", fun = mean, width = 0.6, color = "black") +
  geom_jitter(width = 0.15, size = 4, alpha = 0.7, color = "black") +
  stat_summary(fun.data = mean_se, geom = "errorbar", 
               width = 0.2, linewidth = 1, color = "black") +
  scale_fill_manual(values = survival_colors) +
  theme_pubr(base_size = 16) +
  coord_cartesian(ylim = y_lim) +
  labs(y = "Inverse Simpson Index",
    x = NULL,
    fill = "Survival Status") +
   stat_compare_means(
    aes(group = survival), 
    method = "wilcox.test", 
    label = "p.format", 
    label.y = max(y_lim) * 0.95,
    size = 5) +
    theme(axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 14),
    axis.title.y = element_text(size = 16),
    legend.position = "right",
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    aspect.ratio = 1.5)
p2

# Save as a pdf
pdf("7.1_Post-transplant_WT_only_samples.pdf", width = 6, height = 8)
p2
dev.off()


# # Peter's plots 250214
# p1 <- joined_tibble %>% group_by(patient_id, survival) %>%
#   dplyr::summarize(inv.simpson = mean(inv.simpson)) %>% # take the average for P01
#   ggplot(aes(x = survival, y = inv.simpson)) + 
#   geom_bar(stat = "summary", fun=mean, aes(fill = survival), color = "black") +
#   scale_fill_manual(values=c("skyblue1", "salmon"))+
#   geom_jitter(color = "#00000080", size = 4) +
#   stat_summary(fun.data=mean_se, geom="errorbar", width=.5, linewidth=0.5) +
#   stat_compare_means(aes(group = survival), method = "t.test", label = "p.format", size = 8,
#                      label.y = max(joined_tibble$inv.simpson)*0.95) +
#   ylab("Inverse Simpson Index") +
#   theme_pubr() +
#   theme(strip.text = element_text(size = 20 , color = "black", face="bold"),
#         aspect.ratio = 2,
#         axis.text.x = element_text(angle = 45, vjust= 1, hjust = 1, size = 24, color = "black"),
#         axis.title.x = element_blank(), 
#         axis.text.y = element_text(size = 24),
#         axis.title.y = element_text(size = 24, color = "black"),
#         plot.title =  element_text(size = 24,color = "black", face = "bold"),
#         legend.key.size = unit(10,"mm"),
#         legend.title = element_text(size = 20),
#         legend.position = "right",
#         legend.text = element_text(size = 20))
# 
# ggsave("6.1_AllPatients.pdf", width = 6, height = 8)
# 
# joined_tibble %>% group_by(patient_id, survival) %>% filter(patient_id %in% mt_patients) %>%
#   dplyr::summarize(inv.simpson = mean(inv.simpson)) %>% # take the average for P01
#   ggplot(aes(x = survival, y = inv.simpson)) + 
#   geom_bar(stat = "summary", fun=mean, aes(fill = survival), color = "black") +
#   scale_fill_manual(values=c("skyblue1", "salmon"))+
#   geom_jitter(color = "#00000080", size = 4) +
#   stat_summary(fun.data=mean_se, geom="errorbar", width=.5, linewidth=0.5) +
#   stat_compare_means(aes(group = survival), method = "t.test", label = "p.format", size = 8,
#                       label.y = max(joined_tibble$inv.simpson)*0.9) +
#   ylab("Inverse Simpson Index") +
#   theme_pubr() +
#   theme(strip.text = element_text(size = 20 , color = "black", face="bold"),
#         aspect.ratio = 2,
#         axis.text.x = element_text(angle = 45, vjust= 1, hjust = 1, size = 24, color = "black"),
#         axis.title.x = element_blank(), 
#         axis.text.y = element_text(size = 24),
#         axis.title.y = element_text(size = 24, color = "black"),
#         plot.title =  element_text(size = 24,color = "black", face = "bold"),
#         legend.key.size = unit(10,"mm"),
#         legend.title = element_text(size = 20),
#         legend.position = "right",
#         legend.text = element_text(size = 20))
# ggsave("6.1_TP53mut-only.pdf", width = 6, height = 8)


