# Nurefsan and Peter, 240117-250607
# Plot TCR diversity over time

# Load the libraries
library(tidyverse)
library(Seurat)
library(scRepertoire)
library(cowplot)

# Empty environment
rm(list=ls())

# For Nurefsan:
setwd("/Users/dz855/Dropbox (Partners HealthCare)/ImmuneEscapeTP53/")
# For Peter:
#setwd("~/DropboxMGB/Projects/ImmuneEscapeTP53/TP53_ImmuneEscape/7_TCR_Diversity")

# Favorite function
cutf <- function(x, f=1, d="/") sapply(strsplit(x, d), function(i) paste(i[f], collapse=d))

# Load final Seurat object including TCR calls
seu <- readRDS("../AuxiliaryFiles/250528_Seurat_complete.rds")


######## PREPARE DATA ########

# Keep only annotated T cell clusters
TCAT_cells <- colnames(seu)[!is.na(seu$TCAT_Multinomial_Label)]
seu_T <- subset(seu, cells = TCAT_cells)

# Turn metadata to a tibble and keep only needed variables
metadata_tib <- as_tibble(seu_T@meta.data, rownames = "cell")
metasubset_tib <- metadata_tib %>% select(cell, orig.ident, cohort, patient_id, timepoint, sample_status, TP53_status, CTstrict) 

# Add a new ID column from which to calculate diversity
metasubset_tib$status_timepoint <- paste0(metasubset_tib$sample_status,
  case_when(metasubset_tib$timepoint == 0 ~ "",
            between(metasubset_tib$timepoint, 3, 6) ~ "_3-6M",
            metasubset_tib$timepoint > 6 ~ "_>6M",
            .default = paste0("_", as.character(metasubset_tib$timepoint), "M")))
# Factorize & check
metasubset_tib$status_timepoint <- factor(metasubset_tib$status_timepoint,
    levels = c("pre-transplant", "remission_1M", "remission_1.5M", "remission_3-6M", "remission_>6M", "relapse_3-6M", "relapse_>6M"))
metasubset_tib %>% group_by(sample_status, timepoint, status_timepoint) %>% count()

# Add patient ID, status, and timepoint together as group to calculate diversity by
metasubset_tib$group <- paste0(metasubset_tib$patient_id, "_", metasubset_tib$status_timepoint)

# Select TP53 mutated patients, exclude early relapse cohort
metasubset2_tib <- metasubset_tib %>% filter(TP53_status == "MUT",
  ! patient_id %in% c("P30", "P31", "P32", "P33"))

# Subset to samples with a reasonable number of cells
metasubset2_tib %>% group_by(group) %>% count()
metasubset2_tib %>% group_by(group) %>% count() %>%
  ggplot(aes(x = group, y = n)) +
  geom_point() +
  geom_hline(yintercept = 500, color = "red") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank())

# Filter for samples with >500 cells
metasubset3_tib <- metasubset2_tib %>% group_by(group) %>% filter(n() >= 500)

# Split into list
metasubset_ls <- split(metasubset3_tib, f = metasubset3_tib$group)


######## DIVERSITY FUNCTIONS ########

# This section calculates diversity indices differently from the scRepertoire clonal diversity function.
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
  # df_list <- metasubset_ls
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
  mat$group = names(df_list)
                    
  return(mat)
}


######## CALCULATE AND PLOT DIVERSITY ########

# Calculate the diversity with the function, keep in mind only works with lists 
diversities_df <- compute_diversity(metasubset_ls, "CTstrict", 1000)

# Add information to make annotated plots
metasubset_summary <- metasubset3_tib %>% group_by(group, cohort, patient_id, status_timepoint) %>% summarize
joined_tibble <- diversities_df %>% left_join(metasubset_summary)
  
# Split by cohort
diversities1_df <- joined_tibble %>% subset(subset = cohort == "long-term-remission")
diversities2_df <- joined_tibble %>% subset(subset = cohort == "relapse")

# Parameters for graphics
y_lim <- c(0, max(joined_tibble$inv.simpson))
solid.line.points <- c("pre-transplant", "remission_3-6M")
dashed.line.points <- c("remission_3-6M", "remission_>6M", "relapse_>6M")

# Plot cohort 1
p1 <- ggplot(diversities1_df, aes(x = status_timepoint, y = inv.simpson, group = patient_id)) + 
  geom_point(size = 5, aes(color = patient_id), shape = 1) +
  stat_summary(data = subset(diversities1_df, status_timepoint %in% solid.line.points),
    fun=mean, geom="line", aes(x = status_timepoint, y = inv.simpson, linetype = "solid line", color = patient_id)) +
  stat_summary(data = subset(diversities1_df, status_timepoint %in% dashed.line.points),
    fun=mean, geom="line", aes(x = status_timepoint, y = inv.simpson, linetype = "dashed line", color = patient_id))+
  guides(linetype = "none") +
  scale_color_brewer(palette = "Dark2")+
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
p2 <- ggplot(diversities2_df, aes(x = status_timepoint, y = inv.simpson, group = patient_id)) +
  geom_point(size = 5, aes(color = patient_id), shape = 1) +
  stat_summary(data = subset(diversities2_df, status_timepoint %in% solid.line.points),
    fun = mean, geom = "line", aes(x = status_timepoint, y = inv.simpson, linetype = "solid line", color = patient_id)) +
  stat_summary(data=subset(diversities2_df, status_timepoint %in% dashed.line.points),
    fun = mean, geom = "line", aes(x = status_timepoint, y = inv.simpson, linetype = "dashed line", color = patient_id))+
  guides(linetype = "none") +
  scale_color_brewer(palette = "Set1")+
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

p <- plot_grid(p1,p2)
save_plot("7.4_Longitudinal_diversity.pdf", p, ncol = 2, base_asp = 2, base_height = 4, base_width = 6)
dev.off()



