# name: Nurefsan Sariipek, Peter van Galen
# date: "July 5th, 2023"

library(scRepertoire)
library(Seurat)
library(randomcoloR)
library(RColorBrewer)
library(ggpubr)
library(tidyverse)
library(janitor)
library(ggrepel)
library(cowplot)
library(TSP)
library(ggnewscale)
library(ggrastr)

# Empty environment
rm(list=ls())

# Load the csv files only from cohort 1 and 2 from 2-6mo remission period
# For Nurefsan:
my_wd <- "/Users/dz855/Dropbox (Partners HealthCare)/ImmuneEscapeTP53/"
# For Peter:
#my_wd <- "~/DropboxMGB/Projects/ImmuneEscapeTP53/"

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

# Make a list 
contig_list <- list(P01_1Rem, P01_1RemT, P01_2Rem, P02_1Rem, P04_1Rem, P04_1RemT,
                    P05_1Rem, P06_1Rem, P07_1Rem, P07_1RemT, P08_1Rem, P08_1RemT)

# Combine all the samples
combined <- combineTCR(contig_list,
                       samples = c("P01_1Rem", "P01_1RemT", "P01_2Rem", "P02_1Rem", "P04_1Rem", "P04_1RemT",
                                   "P05_1Rem","P06_1Rem", "P07_1Rem", "P07_1RemT", "P08_1Rem", "P08_1RemT"))

# Add variables
combined <- addVariable(combined, name = "cohort",
                        variables = c("cohort1","cohort1","cohort1","cohort1","cohort1","cohort1",
                                      "cohort2","cohort2","cohort2","cohort2","cohort2","cohort2"))

combined <- addVariable(combined, name = "ptnumber",
                        variables = c("P01_1Rem", "P01_1Rem", "P01_2Rem", "P02_1Rem", "P04_1Rem", "P04_1Rem", "P05_1Rem","P06_1Rem", "P07_1Rem", "P07_1Rem", "P08_1Rem", "P08_1Rem"))


# Load the Seurat object subsetted for T cells
Tcells <- readRDS(paste0(my_wd, "AnalysisNurefsan/RDS files/Tcellsfinal.rds"))
## UMAP dimensions are lost in the previous one check this one
#Tcells <- readRDS(paste0(my_wd, "Trash/old RDS files/Tcellsubset.rds"))

# Keep only annotated T cell clusters (remove NK cells)
Tcells <- subset(x = Tcells, subset = seurat_clusters %in% c(0,1,2,3,4,5,6,7,9,10,11,12,14)) 
# Check the metadata 
#View(Tcells@meta.data)
# Check the T cell types
table(Idents(Tcells))

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

# Select only CD8 
Tcells <- subset(x = Tcells, subset = celltype %in% c("CD8 Effector","CD8 Memory","CD8 Naïve","CD8 Terminally Exhausted"))

# Select only CD4
Tcells <- subset(x = Tcells, subset = celltype %in% c("CD4 Memory","CD4 Naïve","Treg"))

# Plot UMAPs --------------------------------------------------------------------------------------

# Load all data ---------------------------------
AllCells <- readRDS("../../../AnalysisNurefsan/GEX data/seu_diet_merged.rds")

# Define color scheme ---------------------------

# From https://stackoverflow.com/questions/15282580/how-to-generate-a-number-of-most-distinctive-colors-in-r
qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pdf("230724_Colorwheel.pdf")
pie(rep(1,74), col=col_vector)
dev.off()

mycol_tib <- tribble(~celltype, ~color,
                 "HSPCs", col_vector[8],
                 "Early Erythroids", col_vector[36],
                 "Mid Erythroids", col_vector[21],
                 "Late Erythroids", col_vector[22],
                 "Pro Monocytes", col_vector[73],
                 "Monocytes", col_vector[74],
                 "Non Classical Monocytes", col_vector[3],
                 "cDC", col_vector[69],
                 "pDC", col_vector[55],
                 "Pro B cells", col_vector[40],
                 "Pre B cells", col_vector[25],
                 "B cells", col_vector[57],
                 "Plasma Cells", col_vector[5],
                 "CD4 Naïve", col_vector[23],
                 "CD4 Memory", col_vector[28],
                 "Treg", col_vector[11],
                 "CD8 Naïve", col_vector[67],
                 "CD8 Memory", col_vector[63],
                 "CD8 Effector", col_vector[58],
                 "CD8 Terminally Exhausted", col_vector[60],
                 "γδ T lymphocytes", col_vector[1],
                 "NK T cells", col_vector[46],
                 "CD56 Dim NK cells", col_vector[52],
                 "CD56 Bright NK cells", col_vector[61],
                 "Blasts", col_vector[45])
mycol <- mycol_tib$color
names(mycol) <- mycol_tib$celltype

# Plot UMAP with all cells ----------------------
p1 <- as_tibble(AllCells@meta.data) %>%
  mutate(UMAP1 = AllCells@reductions$umap@cell.embeddings[,1],
         UMAP2 = AllCells@reductions$umap@cell.embeddings[,2]) %>%
  filter(! celltype %in% c("Unidentified", "Doublets")) %>%
  mutate(celltype = gsub("CD8 Central Memory", "CD8 Memory", gsub("CD8 Effector Memory", "CD8 Effector", celltype))) %>%
  mutate(celltype = factor(celltype, levels =
    c("HSPCs", "Early Erythroids", "Mid Erythroids", "Late Erythroids", 
    "Pro Monocytes", "Monocytes", "Non Classical Monocytes", "cDC", "pDC",
    "Pro B cells", "Pre B cells", "B cells", "Plasma Cells",
    "CD4 Naïve", "CD4 Memory", "Treg",
    "CD8 Naïve", "CD8 Memory", "CD8 Effector", "CD8 Terminally Exhausted",
    "γδ T lymphocytes",
    "NK T cells", "CD56 Dim NK cells", "CD56 Bright NK cells", 
    "Blasts"))) %>%
  ggplot(aes(x = -UMAP1, y = UMAP2, color = celltype)) +
  geom_point_rast(size = 0.1) +
  scale_color_manual(values = mycol) +
  theme_pubr() +
  guides(color = guide_legend(override.aes = list(size = 2.5))) +
  theme(aspect.ratio = 1, legend.position = "right")

pdf("230724_UMAP_all.pdf", width = 8, height = 5)
print(p1)
dev.off()

# Plot UMAP with T cells ------------------------
p2 <- as_tibble(Tcells@meta.data) %>%
  mutate(UMAP1 = Tcells@reductions$umap@cell.embeddings[,1],
         UMAP2 = Tcells@reductions$umap@cell.embeddings[,2]) %>%
  mutate(celltype = factor(celltype, levels =
                             c("CD4 Naïve", "CD4 Memory", "Treg",
                               "CD8 Naïve", "CD8 Memory", "CD8 Effector", "CD8 Terminally Exhausted",
                               "γδ T lymphocytes", 
                               "NK T cells"))) %>%
  ggplot(aes(x = UMAP1, y = UMAP2, color = celltype)) +
  geom_point_rast(size = 0.1) +
  scale_color_manual(values = mycol) +
  theme_pubr() +
  guides(color = guide_legend(override.aes = list(size = 2.5))) +
  theme(aspect.ratio = 1, legend.position = "right")


pdf("230724_UMAP_T.pdf", width = 8, height = 5)
p2
dev.off()

# Other visualization ---------------------------

# Add color palette (not subsequently used)
colorblind_vector <- colorRampPalette(rev(c("#0D0887FF", "#47039FFF", 
                                            "#7301A8FF", "#9C179EFF", "#BD3786FF", "#D8576BFF",
                                            "#ED7953FF","#FA9E3BFF", "#FDC926FF", "#F0F921FF")))

# Combine TCR and sc-RNAseq data
Tcells_combined <- combineExpression(combined, Tcells, 
                                     cloneCall = "strict",
                                     group.by = "ptnumber",
                                     proportion = FALSE,
                                     filterNA = T,
                                     cloneTypes = c(Rare = 1e-04, Small = 0.001, Medium = 0.01, Large = 0.1, Hyperexpanded =1))

# save this combined object to use in other purposes.
saveRDS(Tcells_combined, file = "combined.RDS")

# Calculate the frequency, setting group by to sample which is combined samples(T-cell enriched and MNC)
clonalDiversity(Tcells_combined,
                cloneCall = "strict",
                group.by = "sample",
                metrics = c("inv.simpson","gini.simpson"),
                exportTable = T)

# Adding Souporcell information

# Load the metadata that contains souporcell information
combined_df <- read_csv(paste0(my_wd, "/AnalysisNurefsan/Souporcell/output/cohort1-2_souporcell.csv"))

# Wrangle the metadata to 
combined_df$cell = gsub("_.*","", combined_df$cell)
combined_df$id = gsub("\\.","_",combined_df$id )
combined_df$cell= paste0(combined_df$id, "_",combined_df$cell)

# Left join the 2 metadata
Tcells_combined_tib <- as_tibble(Tcells_combined@meta.data, rownames = "cell")
newdf <- Tcells_combined_tib %>% 
  left_join(combined_df, by ="cell") %>% 
  drop_na()

# Add assignment calls to the Seurat metadata
Tcells_combined <- AddMetaData(Tcells_combined, data.frame(select(newdf, cell, assignment), row.names = "cell"))
# Check the numbers
Tcells_combined$assignment %>% tabyl

# Look into different CD8 T cells
Tcells_combined_cd8_Effector <- subset(x = Tcells_combined, subset = celltype == "CD8 Effector")
Tcells_combined_cd8_Memory <- subset(x = Tcells_combined, subset = celltype == "CD8 Memory")
Tcells_combined_cd8_naive <- subset(x = Tcells_combined, subset = celltype == "CD8 Naïve")
Tcells_combined_cd8_exhausted <- subset(x = Tcells_combined, subset = celltype =="CD8 Terminally Exhausted")


# Subset only donor cells 
Tcells_combined_donor <- subset(x = Tcells_combined, subset = assignment == "donor")
Tcells_combined_donor_cd8 <- subset(x = Tcells_combined_donor, subset = celltype %in% c("CD8 Effector","CD8 Memory","CD8 Naïve","CD8 Terminally Exhausted"))
Tcells_combined_donor_cd4 <- subset(x = Tcells_combined_donor, subset = celltype %in% c("CD4 Memory","CD4 Naïve","Treg"))

# Subset only host cells 
Tcells_combined_host <- subset(x = Tcells_combined, subset = assignment == "host")
Tcells_combined_host_cd8 <- subset(x = Tcells_combined_host, subset = celltype %in% c("CD8 Effector","CD8 Memory","CD8 Naïve","CD8 Terminally Exhausted"))
Tcells_combined_host_cd4 <- subset(x = Tcells_combined_host, subset = celltype %in% c("CD4 Memory","CD4 Naïve","Treg"))

# Calculate the inverse simpson index for each 

clonalDiversity(Tcells_combined_cd8_exhausted,
                cloneCall = "strict",
                group.by = "sample",
                metrics = c("inv.simpson","gini.simpson"),
                exportTable = T)

clonalDiversity(Tcells_combined_cd8,
                cloneCall = "strict",
                group.by = "sample",
                metrics = c("inv.simpson","gini.simpson"),
                exportTable = T)

clonalDiversity(Tcells_combined_cd8,
                cloneCall = "strict",
                group.by = "id",
                x.axis = "celltype",
                metrics = c("inv.simpson","gini.simpson"))

# Recalculate Frequency ---------------------------------------------------------------------------

#As of 231226 this section is unneccessary since we have corrected the way that screpertoire calculates the frequency and now both our calculation and screpertoires match up.

# The combineExpression function from scRepertoire added a Frequency column to the object, which is
# paramount. Here, I will recalculate the Frequency and compare it to the scRepertoire value.

# First, create a tibble of the metadata and add UMAP coordinates
Tcells_combined_tib <- as_tibble(Tcells_combined@meta.data, rownames = "cell")
Tcells_combined_tib$UMAP1 <- Tcells_combined@reductions$umap@cell.embeddings[,1]
Tcells_combined_tib$UMAP2 <- Tcells_combined@reductions$umap@cell.embeddings[,2]

# Add one column to identify the patient/timepoint (pool MNC and CD3-enriched libraries) and one
# column to group by clonotype therein
Tcells_combined_tib <- Tcells_combined_tib %>%
  mutate(pt_timepoint = gsub("RemT", "Rem", id)) %>%
  mutate(pt_timepoint_ct = paste0(pt_timepoint, "_", CTstrict))

# Now calculate the Frequency / clone size
Tcells_combined_tib <- Tcells_combined_tib %>% group_by(pt_timepoint_ct) %>% add_count() %>%
  ungroup() %>% arrange(pt_timepoint_ct) %>% arrange(pt_timepoint, n)
# Wiew(Tcells_combined_tib)

# Visualizations the clonotype size
Tcells_combined_tib %>% mutate(rown = row_number()) %>%
  ggplot(aes(x = rown, y = n, color = pt_timepoint)) +
  geom_point() +
  theme_pubr()

Tcells_combined_tib %>% mutate(rown = row_number()) %>%
  ggplot(aes(x = rown, y = n, color = cohort)) +
  geom_point() +
  theme_pubr()

# Compare scRepertoire: Frequency, my code: n
Tcells_combined_tib %>%
  ggplot(aes(x = clonalFrequency, y = n, color = pt_timepoint)) +
  geom_point() +
  theme_pubr() +
  theme(aspect.ratio = 1)
# That's pretty different! I think my code is more correct, because it pools MNC and T cell
# enriched samples, and because it calculates clonotype size within pt_timepoint groups.
# However, it is still confounded by the number of cells.

# The cell numbers are different per cohort and per patient, which confounds subsequent analysis
# After considering subsetting to 10,000 T cells per cohort ...
#Tcells_subset_tib <- Tcells_combined_tib %>% group_by(cohort) %>% sample_n(size = 10000) %>% ungroup
# ...or 220 CD8+ T cells per pt_timepoint...
#Tcells_subset_tib <- Tcells_combined_tib %>% filter(grepl("CD8", celltype)) %>%
#  group_by(pt_timepoint) %>% sample_n(size = 220) %>% ungroup
# ...I decided this is the most fair solution:
Tcells_combined_tib %>% tabyl(pt_timepoint)
Tcells_subset_tib <- Tcells_combined_tib %>% group_by(pt_timepoint) %>% slice_sample(n = 513) %>%
  select(-n) %>% group_by(pt_timepoint_ct) %>% add_count() %>%
  ungroup() %>% arrange(pt_timepoint_ct) %>% arrange(pt_timepoint, n)
# The sampling here matters quite a bit - we should do the average of 1000 samplings or sth like that

# Visualize subsampled data
Tcells_subset_tib %>%
  mutate(rown = row_number()) %>%
  ggplot(aes(x = rown, y = n, color = pt_timepoint)) +
  geom_point() +
  theme_pubr()


# Cumulative plot ---------------------------------------------------------------------------------

# Summarize by clonotype, then exclude clones of size 1
Tcells_subset_summary_tib <- Tcells_subset_tib %>% group_by(cohort, pt_timepoint, CTstrict) %>% summarize(n = n())
Tcells_subset_summary_tib <- Tcells_subset_summary_tib %>% filter(n >= 2)

# Add normalized frequency and normalized clone index
Tcells_subset_summary_tib <- Tcells_subset_summary_tib %>% group_by(pt_timepoint) %>%
  arrange(pt_timepoint, -n) %>% mutate(n_norm = n/sum(n)) %>% mutate(n_cumsum = cumsum(n_norm))
Tcells_subset_summary_tib <- Tcells_subset_summary_tib %>%
  mutate(clone_index = row_number()) %>%
  mutate(clone_index_norm = clone_index / max(clone_index))

# Visualize
p1 <- Tcells_subset_summary_tib %>%
  ggplot(aes(x = clone_index_norm, y = n_cumsum, color = cohort, group = pt_timepoint)) +
  geom_line(size = 1) +
  theme_pubr() +
  xlab("Clones ranked by size") +
  ylab("Cumulative proportion") +
  scale_color_manual(values = c("#00FF40CC", "#FF0800CC")) +
  theme(aspect.ratio = 1,
        axis.ticks.x = element_blank(), axis.text.x = element_blank())

# Color by patient
p2 <- Tcells_subset_summary_tib %>%
  ggplot(aes(x = clone_index_norm, y = n_cumsum, color = pt_timepoint, group = pt_timepoint)) +
  geom_line(size = 1) +
  theme_pubr() +
  xlab("Clones ranked by size") +
  ylab("Cumulative proportion") +
  theme(aspect.ratio = 1,
        axis.ticks.x = element_blank(), axis.text.x = element_blank())

# Save
pdf("230723_Cumulative_proportion.pdf", width = 12, height = 5)
plot_grid(p1, p2)
dev.off()


# Connectivity plot -------------------------------------------------------------------------------

# Function to find the shortest path between connected points
sort_func <- function(umap1 = "UMAP1", umap2 = "UMAP2") {
  umap_coord <- data.frame(umap1 = umap1, umap2 = umap2)
  dist_matrix <- dist(umap_coord)
  tour <- solve_TSP(TSP(dist_matrix))
  as.vector(tour)
}

# Split into a list of tibbles, then calculate the shortest path using the Traveling Salesperson Problem algorithm
# The first line is required to ensure groups stay in the same order
Tcells_subset_tib$pt_timepoint_ct <- factor(Tcells_subset_tib$pt_timepoint_ct,
                                            levels = unique(Tcells_subset_tib$pt_timepoint_ct))
clonotypes_ls <- split(Tcells_subset_tib, Tcells_subset_tib$pt_timepoint_ct)
# This speeds up subsequent steps
clonotypes_ls <- lapply(clonotypes_ls, data.frame)
# Calculate shortest path
tours <- lapply(clonotypes_ls, sort_func)
# Then order by the path and recreate a dataframe
clonotypes_order_ls <- lapply(seq_along(clonotypes_ls), function(x) clonotypes_ls[[x]][tours[[x]],]  )
clonotypes_order_tib <- as_tibble(do.call(rbind, clonotypes_order_ls))
# Check it all stayed in order
plot(Tcells_subset_tib$n)
points(clonotypes_order_tib$n, col = "red", pch = ".")

# Use same UMAP coordinates and size scale for both plots
ulimits <- c(min(Tcells_combined_tib$UMAP1), max(Tcells_combined_tib$UMAP1),
             min(Tcells_combined_tib$UMAP2), max(Tcells_combined_tib$UMAP2))
slimits <- c(0, max(clonotypes_order_tib$n))

# Wrangle cohort 1
cohort1_tib <- clonotypes_order_tib %>% filter(cohort == "cohort1") %>%
  mutate(CTstrict_numeric = as.numeric(as.factor(CTstrict))) %>%
  mutate(clone_id = paste0(pt_timepoint, ", Clonotype ", CTstrict_numeric, ", ", n, " cells"))
top_clones1 <- cohort1_tib %>% arrange(-n) %>%
  filter(clone_id %in% unique(clone_id)[1:5])
# SKIP FOR CONNECTIVITY-3
cohort1_tib <- cohort1_tib %>% mutate(clone_id = ifelse(clone_id %in% top_clones1$clone_id, yes = clone_id, no = NA))

# Wrangle cohort 2
cohort2_tib <- clonotypes_order_tib %>% filter(cohort == "cohort2") %>%
  mutate(CTstrict_numeric = as.numeric(as.factor(CTstrict))) %>%
  mutate(clone_id = paste0(pt_timepoint, ", Clonotype ", CTstrict_numeric, ", ", n, " cells")) 
top_clones2 <- cohort2_tib %>% arrange(-n) %>%
  filter(clone_id %in% unique(clone_id)[1:5])
# SKIP FOR CONNECTIVITY-3
cohort2_tib <- cohort2_tib %>% mutate(clone_id = ifelse(clone_id %in% top_clones2$clone_id, yes = clone_id, no = NA))


# CONNECTIVITY-1 (original version)

p1 <- ggplot() +
  geom_path(data = top_clones1,
            aes(x = UMAP1, y = UMAP2, group = clone_id, color = clone_id),
            linewidth = 0.5) +
  new_scale_color() +
  geom_point(data = cohort1_tib,
             aes(x = UMAP1, y = UMAP2, group = pt_timepoint_ct, color = n),
             size = 0.5) +
  scale_color_continuous(limits = slimits) +
  coord_cartesian(xlim = ulimits[1:2], ylim = ulimits[3:4]) +
  ggtitle("Cohort 1 (remission)") +
  labs(color = "Clonotype size") +
  theme_pubr() +
  theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5), legend.position = "right")

p2 <- ggplot() +
  geom_path(data = top_clones2,
            aes(x = UMAP1, y = UMAP2, group = clone_id, color = clone_id),
            linewidth = 0.5) +
  new_scale_color() +
  geom_point(data = cohort2_tib,
             aes(x = UMAP1, y = UMAP2, group = pt_timepoint_ct, color = n),
             size = 0.5) +
  scale_color_continuous(limits = slimits) +
  coord_cartesian(xlim = ulimits[1:2], ylim = ulimits[3:4]) +
  ggtitle("Cohort 2 (relapse)") +
  labs(color = "Clonotype size") +
  theme_pubr() +
  theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5), legend.position = "right")

pdf("230723_Connectivity.pdf", width = 12, height = 4)
plot_grid(p1, p2)
dev.off()


# CONNECTIVITY-2 (connecting lines only)

p1 <- cohort1_tib %>% mutate(clone_id = factor(clone_id, levels = unique(arrange(cohort1_tib, -n)$clone_id))) %>%
  ggplot(aes(x = UMAP1, y = UMAP2, group = clone_id, color = clone_id)) +
  geom_point(size = 1.5) +
  #geom_path(data = top_clones1, linewidth = 1) +
  #scale_color_manual(values = col_vector[c(5, 12, 60, 49, 13)], na.value = "#0000004D") +
  coord_cartesian(xlim = ulimits[1:2], ylim = ulimits[3:4]) +
  ggtitle("Cohort 1 (remission)") +
  labs(color = "Cohort 1 five largest clonotypes") +
  theme_pubr() +
  theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5), legend.position = "right")

p2 <- cohort2_tib %>% mutate(clone_id = factor(clone_id, levels = unique(arrange(cohort2_tib, -n)$clone_id))) %>%
  ggplot(aes(x = UMAP1, y = UMAP2, group = clone_id, color = clone_id)) +
  geom_point(size = 1.5) +
  #geom_path(data = top_clones2, linewidth = 1) +
  #scale_color_manual(values = col_vector[c(53, 26, 51, 48, 47)], na.value = "#0000004D") +
  coord_cartesian(xlim = ulimits[1:2], ylim = ulimits[3:4]) +
  ggtitle("Cohort 2 (relapse)") +
  labs(color = "Cohort 2 five largest clonotypes") +
  theme_pubr() +
  theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5), legend.position = "right")

pdf("230727_Connectivity.pdf", width = 12, height = 4)
plot_grid(p1, p2)
dev.off()



# CONNECTIVITY-3 (color/connect clonotypes with > 5 cells)

p1 <- cohort1_tib %>%
  ggplot(aes(x = UMAP1, y = UMAP2, group = clone_id, color = clone_id)) +
  geom_point(data = filter(cohort1_tib, n <= 5), size = 1.5, color = "black") +
  geom_path(data = filter(cohort1_tib, n > 5), linewidth = 0.5) +
  geom_point(data = filter(cohort1_tib, n > 5), size = 1.5) +
  scale_color_manual(values = sample(col_vector, 2000, replace = T)) +
  coord_cartesian(xlim = ulimits[1:2], ylim = ulimits[3:4]) +
  ggtitle("Cohort 1 (remission)") +
  theme_pubr() +
  theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5), legend.position = "none")

p2 <- cohort2_tib %>%
  ggplot(aes(x = UMAP1, y = UMAP2, group = clone_id, color = clone_id)) +
  geom_point(data = filter(cohort2_tib, n <= 5), size = 1.5, color = "black") +
  geom_path(data = filter(cohort2_tib, n > 5), linewidth = 0.5) +
  geom_point(data = filter(cohort2_tib, n > 5), size = 1.5) +
  scale_color_manual(values = sample(col_vector, 2000, replace = T)) +
  coord_cartesian(xlim = ulimits[1:2], ylim = ulimits[3:4]) +
  ggtitle("Cohort 2 (relapse)") +
  theme_pubr() +
  theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5), legend.position = "none")

pdf("230728_Connectivity.pdf", width = 12, height = 4)
plot_grid(p1, p2)
dev.off()



# CONNECTIVITY-4 (color top 500 cells)
cohort1_tib <- cohort1_tib %>% arrange(-n) %>% mutate(top = c(rep("Yes", 172), rep("No", nrow(cohort1_tib)-172)))
cohort2_tib <- cohort2_tib %>% arrange(-n) %>% mutate(top = c(rep("Yes", 172), rep("No", nrow(cohort2_tib)-172)))

# How many clones are in the top 172 cells?
cohort1_tib %>% slice_head(n = 172) %>% .$clone_id %>% unique %>% length
cohort2_tib %>% slice_head(n = 172) %>% .$clone_id %>% unique %>% length


p1 <- cohort1_tib %>%
  ggplot(aes(x = UMAP1, y = UMAP2, group = clone_id, color = clone_id)) +
  geom_point(data = filter(cohort1_tib, top == "No"), size = 1.5, color = "black") +
  geom_path(data = filter(cohort1_tib, top == "Yes"), linewidth = 0.5) +
  geom_point(data = filter(cohort1_tib, top == "Yes"), size = 1.5) +
  scale_color_manual(values = sample(col_vector)) +
  coord_cartesian(xlim = ulimits[1:2], ylim = ulimits[3:4]) +
  ggtitle("Cohort 1 (remission)") +
  theme_pubr() +
  theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5), legend.position = "none")

p2 <- cohort2_tib %>%
  ggplot(aes(x = UMAP1, y = UMAP2, group = clone_id, color = clone_id)) +
  geom_point(data = filter(cohort2_tib, top == "No"), size = 1.5, color = "black") +
  geom_path(data = filter(cohort2_tib, top == "Yes"), linewidth = 0.5) +
  geom_point(data = filter(cohort2_tib, top == "Yes"), size = 1.5) +
  scale_color_manual(values = sample(col_vector)) +
  coord_cartesian(xlim = ulimits[1:2], ylim = ulimits[3:4]) +
  ggtitle("Cohort 2 (relapse)") +
  theme_pubr() +
  theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5), legend.position = "none")

pdf("230728_Connectivity-2.pdf", width = 12, height = 4)
plot_grid(p1, p2)
dev.off()



# THE FOLLOWING NEEDS TO BE REASSESSED

# Correlate gene expression to clonotype size -----------------------------------------------------

# Exctract gene expression data
expr_mat <- as.matrix(LayerData(Tcells_combined, layer = "data"))

Tcells_combined_tib

expr_mat[1:3,1:3]
dim(expr_mat)
Tcells_combined_tib$n
# For each gene, determine correlation and p-value to clonotype size (Frequency)
gene_freq_cor <- apply(expr_mat, 1, function(x) cor(x, meta_T$Frequency))
gene_freq_p <- apply(expr_mat, 1, function(x) cor.test(x, meta_T$Frequency)$p.value)

# Make a tibble with these values
gene_freq_cor_tib <- tibble(gene = names(gene_freq_cor), correlation = gene_freq_cor, pvalue = gene_freq_p)
gene_freq_cor_tib <- na.omit(gene_freq_cor_tib) %>% arrange(-correlation)

# Select genes to highlight: from Penter et al, 2021
penter_exh <- c("PDCD1", "HAVCR2", "LAG3", "TIGIT", "KLRG1", "CD38", "TOX", "TBX21", "EOMES", "CD244", "TNFRSF9", "GZMB", "PRF1")
# Add cherrypicked genes from the top of the list (think about addinge more, including negative?)
print(gene_freq_cor_tib, n = 50)
highlight_genes <- c("NKG7", "CCL5", "CD8A", "CD8B", "CD274", "GZMA", penter_exh)

# Visualize
gene_freq_cor_tib %>%
  ggplot(aes(x = correlation, y = -log10(pvalue), label = gene)) +
  geom_point() +
  geom_text_repel(data = subset(gene_freq_cor_tib, gene %in% highlight_genes), min.segment.length = 0.1, color = "red") +
  theme_pubr() +
  theme(aspect.ratio = 1)


