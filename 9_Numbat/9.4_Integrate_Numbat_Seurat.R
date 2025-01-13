# Nurefsan Sariipek, 241219
# Combine Numbat results with gen-ex data

# Load Libraries
library(readr)
library(numbat)
library(Seurat)
library(tidyverse)
library(RColorBrewer)
library(ggsci)

# Empty environment
rm(list=ls())

# Set working directory
# For Nurefsan
setwd("/Users/dz855/Dropbox (Partners HealthCare)/ImmuneEscapeTP53/TP53_ImmuneEscape/9_Numbat")
# For Peter
#setwd("~/DropboxMGB/Projects/ImmuneEscapeTP53/TP53_ImmuneEscape/9_Numbat")

# Load the saved Seurat objects
seu_diet_merged <- readRDS("../../RDS files/seu_diet_merged.rds")

# Check which cells are labeled as blast cells 
#blast_percentages <- seu_diet_merged@meta.data %>%
#  group_by(Sample) %>%
#  summarize(
#    total_cells = n(),
#    blast_cells = sum(celltype == "Blasts"),
#    blast_percentage = round((blast_cells / total_cells) * 100, 4))
# View results
#print(blast_percentages)

############################# COMBINE NUMBAT OUTPUT TO SEURAT #############################

# Load the saved dataframe that contains souporcell information for selecting only host counts
# For cohort 1&2
soc_combined_df <- read_csv("../5_Souporcell/results/cohort1-2_souporcell.csv")
# For cohort 3
soc_combined_df <- read_csv("../5_Souporcell/results/cohort3_souporcell.csv")

# Subset souporcell output for mononuclear cell libraries 
# For patient 5
soc_subset <- soc_combined_df %>% filter(orig.ident %in% c("2737_MNC","25809_MNC","9596_MNC"), assignment == "host")

# For patient 8
soc_subset <- soc_combined_df %>% filter(orig.ident %in% c("4618_MNC","6174_MNC","9931_MNC","1953_MNC"), assignment == "host")

# For patient 9
soc_subset <- soc_combined_df %>% filter(orig.ident %in% c("1677_MNC","1732_MNC","1811_MNC"), assignment == "host")

#For patient 10
soc_subset <- soc_combined_df %>% filter(orig.ident %in% c("1195_MNC","1285_MNC","1347_MNC"), assignment == "host")

# For patient 12
soc_subset <- soc_combined_df %>% filter(orig.ident %in% c("9355_MNC","1013_MNC"), assignment == "host")

# Subset Seurat object for desired cells
seu_subset <- subset(seu_diet_merged, cells = soc_subset$cell)

# Modify the barcodes to keep '-1' and remove everything after
colnames(seu_subset) <- sub("-1.*", "-1", colnames(seu_subset))

# Load the TSV file that has clonal information from Numbat for Patient 5
pt5 <- read.table("/Volumes/sariipek/numbat/Pt05/hostpt5/clone_post_2.tsv", row.names = 1, header = T)
pt8 <- read.table("/Volumes/sariipek/numbat/Pt08/hostpt8/clone_post_2.tsv", row.names = 1, header = T)
pt9 <- read.table("/Volumes/sariipek/numbat/Pt09/hostpt9/clone_post_2.tsv", row.names = 1, header = T)
pt10 <- read.table("/Volumes/sariipek/numbat/Pt10/hostpt10_2/clone_post_1.tsv", row.names = 1, header = T)
pt12 <- read.table("/Volumes/sariipek/numbat/Pt12/hostpt12/clone_post_2.tsv", row.names = 1, header = T)
# For Peter
#pt5 <- read.table("/Volumes/broad_vangalenlab/sariipek/numbat/Pt05/hostpt5/clone_post_2.tsv", row.names = 1, header = T)
pt_select  <- pt5[, c("clone_opt", "compartment_opt")]
head(pt_select)

# Check
all(rownames(pt_select) %in% colnames(seu_subset))
all(colnames(seu_subset) %in% rownames(pt_select))
# This is probably unneccessary, but order pt5 in the same way as the Seurat object before merging
pt_order <- pt_select[colnames(seu_subset),]

# Add Numbat results to the Seurat object
seu_subset <- AddMetaData(seu_subset, metadata = pt_select)

# Verify the added metadata
View(seu_subset@meta.data)

# UMAP visualization
DimPlot(seu_subset, reduction = "umap", group.by = "compartment_opt") + theme(aspect.ratio = 1)
DimPlot(seu_subset, reduction = "umap", group.by = "celltype", split.by = "compartment_opt", label = T) + scale_color_igv() + theme(aspect.ratio = 1, legend.position = "none")

#Save this combined seurat object for each patient
saveRDS(seu_subset,"/Users/dz855/Partners HealthCare Dropbox/Nurefsan Sariipek/ImmuneEscapeTP53/TP53_ImmuneEscape/9_Numbat/pt5.RDS")

# End of combination of Numbat and Seurat object

############################# VISUALIZATIONS #############################
# Empty environment
rm(list=ls())

#Load from saved object
seu_subset <- readRDS("/Users/dz855/Partners HealthCare Dropbox/Nurefsan Sariipek/ImmuneEscapeTP53/TP53_ImmuneEscape/9_Numbat/pt9.RDS")

# Make a dataframe of malignant cells according to Numbat
metadata_subset1 <- as_tibble(seu_subset@meta.data, rownames = "cell") 
#%>% subset(compartment_opt == "normal")

# Group and summarize the data
summarized_data <- metadata_subset1 %>%
  group_by(status, celltype, compartment_opt) %>%
  summarise(cells = n(), .groups = "drop")  # Count cells per group

# Re-order the status information for visualization purposes
summarized_data$status <- factor(summarized_data$status, levels = c("pre_transplant", "remission", "relapse"))

# Calculate normalized proportions and total cell numbers
summarized_data <- summarized_data %>%
  group_by(compartment_opt) %>%
  mutate(total_cells = sum(cells)) %>%
  mutate(cells_normalized = cells / total_cells)

  total_cells <- summarized_data %>%
  group_by(compartment_opt) %>%
  summarise(total_cells = sum(cells))

# Create the bar graph to visualize different celltypes across the annotated cells
p1 <- ggplot(summarized_data, aes(x = status, y = cells, fill = celltype)) +
  geom_bar(stat = "identity", position = "stack") +  # Stacked bar chart
  labs(#title = "Cell Counts by Timepoint and Compartment",
    x = "Timepoint",
    y = "Number of Cells",
    fill = "Celltype "
  ) +
  scale_fill_igv() +
  theme_bw()+
  theme(aspect.ratio = 0.75,
                   axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1,face="plain", size=18, color="black"), 
                   axis.text.y = element_text(face="plain", size=16, color="black"),
                   axis.title.y = element_text(size = 18),
                   axis.title.x = element_blank(),
                   plot.title = element_text(size=20, face="plain"),
                   strip.text = element_text(size=12, face="bold"),
                   legend.title=element_text(size=16), 
                   legend.text=element_text(size=16)) +
  theme(panel.grid = element_blank())

p1

pdf("9.4_Pt10_tumorcells.pdf", width = 16, height = 8)
p1
dev.off()

# Normalized tumor cell ratio visualization
p2 <- ggplot(summarized_data, aes(x = compartment_opt, y = cells_normalized, fill = celltype)) +
  geom_bar(stat = "identity", position = "stack") +  # Stacked bar chart
  geom_text(data = total_cells, aes(x = compartment_opt, y = 1.05, label = total_cells),
            inherit.aes = FALSE,  # Do not inherit aesthetics from the main ggplot
            size = 6, color = "black", fontface = "bold") +  
  labs(#title = "Cell Counts by Timepoint and Compartment",
    x = "Timepoint",
    y = "Proportion of Cells (Normalized to 1)",
    fill = "Celltype "
  ) +
  scale_fill_igv() +
  theme_bw()+
  theme(aspect.ratio = 0.75,
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1,face="plain", size=18, color="black"), 
        axis.text.y = element_text(face="plain", size=16, color="black"),
        axis.title.y = element_text(size = 18),
        axis.title.x = element_blank(),
        plot.title = element_text(size=20, face="plain"),
        strip.text = element_text(size=12, face="bold"),
        legend.title=element_text(size=16), 
        legend.text=element_text(size=16),
        panel.grid = element_blank())
p2

pdf("9.4_Pt10_normalcells_normalized.pdf", width = 16, height = 8)
p2
dev.off()

#All cells normalized
p3 <- ggplot(summarized_data, aes(x = compartment_opt, y = cells_normalized, fill = celltype)) +
  geom_bar(stat = "identity", position = "stack") +  # Stacked bar chart
  # Add celltype labels on each stacked segment
  geom_text(aes(label= ifelse(cells_normalized > 0.03, celltype, "")),
            position = position_stack(vjust = 0.5),  # Centered within each stack
            size = 3, 
            color = "black", ) + 
  # Add total cell counts at the top of each bar
  geom_text(data = total_cells, aes(x = compartment_opt, y = 1.05, label = total_cells),
            inherit.aes = FALSE,  # Do not inherit aesthetics from the main ggplot
            size = 6, color = "black", fontface = "bold") +
  labs(#title = "Cell Counts by Timepoint and Compartment",
    x = "Timepoint",
    y = "Proportion of Cells (Normalized to 1)",
    fill = "Celltype") +
  scale_fill_igv() +
  theme_bw()+
  theme(aspect.ratio = 0.75,
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1,face="plain", size=18, color="black"), 
        axis.text.y = element_text(face="plain", size=16, color="black"),
        axis.title.y = element_text(size = 18),
        axis.title.x = element_blank(),
        plot.title = element_text(size=20, face="plain"),
        strip.text = element_text(size=12, face="bold"),
        legend.position = "none",
        panel.grid = element_blank())
p3

pdf("9.4_Pt9_allcells_normalized.pdf", width = 16, height = 8)
p3
dev.off()
