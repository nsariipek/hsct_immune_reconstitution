# Nurefsan Sariipek, 250509
# Load the needed libraries
library(Seurat)
library(dplyr)
library(tidyverse)
library(DESeq2)
library(apeglm)
library(data.table)
library(cowplot)
library(ggplot2)
library(ggpubr)
library(pheatmap)

# Empty environment
rm(list=ls())

# Set working directory
setwd("~/TP53_ImmuneEscape/3_DGE/")

# Load the saved seurat objects
seu <- readRDS("~/250426_Seurat_annotated.rds")

#Load the souporcellresult and add as a metadata
souporcell_df <- read_csv("~/250428_final_dataset.csv")

# Wrangle the df
souporcell_df <- souporcell_df %>%
  mutate(cell = paste(orig.ident, barcode, sep = "_")) 

df_meta <- souporcell_df %>%
  select(cell, origin) %>%
  filter(cell %in% Cells(seu)) %>%            
  distinct(cell, .keep_all = TRUE) %>%        
  column_to_rownames("cell")  

seu <- AddMetaData(seu, metadata = df_meta)

# Normalize
seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)
# Scale the data
seu <- ScaleData(seu)

# Select the group you want to run analysis
ex <- seu %>% subset(cohort=="relapse" & timepoint %in% c("3","5","6") & celltype %in% c("HSC MPP","MEP", "LMPP","Cycling Progenitor", "Early GMP") & sample_status== "remission"& origin %in% c("donor", "recipient"))

# Combine these cell types into one label
ex$celltype_merged <- "Merged_Progenitors"

# Check that how many cells are per pseudo-sample
ex@meta.data %>%
  group_by(sample_id, origin, celltype_merged) %>%
  summarize(n=n()) %>%
  pivot_wider(names_from = "celltype_merged", values_from = "n") %>% View()


# Use all cells in 'ex' (which are now labeled Merged_Progenitors)
seurat_ct <- ex
ct <- "Merged_Progenitors"

# Aggregate to get pseudobulk
bulk_ct <- AggregateExpression(
  seurat_ct,
  return.seurat = TRUE,
  assays = "RNA",
  group.by = c("sample_id", "origin")
)

# Format sample_id
bulk_ct@meta.data$sample_id <- str_replace_all(bulk_ct@meta.data$sample_id, "_", "-")

# Add number of cells per sample
n_cells <- seurat_ct@meta.data %>%
  dplyr::count(sample_id, origin) %>%
  dplyr::rename("n_cells" = "n") %>%
  mutate(orig.ident = paste0(str_replace_all(sample_id, "_", "-"), "_", origin)) %>%
  dplyr::select(-c(origin, sample_id))

meta_bulk_ct <- left_join(bulk_ct@meta.data, n_cells, by = "orig.ident")
rownames(meta_bulk_ct) <- meta_bulk_ct$orig.ident
bulk_ct@meta.data <- meta_bulk_ct

# Exclude unpaired samples
bulk_ct@meta.data$origin <- factor(bulk_ct@meta.data$origin)
if (length(unique(bulk_ct@meta.data$origin)) < 2) stop("Not enough groups to compare.")

# Check the cell counts per group
barplot <- ggplot(bulk_ct@meta.data, aes(x = sample_id, y = n_cells, fill = origin)) +
  geom_bar(stat = "identity", color = "black") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = paste("Cell Counts:", ct), x = "Sample", y = "Number of Cells") +
  geom_text(aes(label = n_cells), vjust = -0.5)
ggsave(paste0("cell_counts_", gsub(" ", "_", ct), ".pdf"), barplot, width = 8, height = 5)

#Run DESeq 
cluster_counts <- FetchData(bulk_ct, layer = "counts", vars = rownames(bulk_ct))

dds <- DESeqDataSetFromMatrix(t(cluster_counts), colData = bulk_ct@meta.data, design = ~ sample_id + origin) # we decided ~ smaple_id +origin is better 
# This line put donor as a reference, but check line 136
dds$origin <- relevel(dds$origin, ref = "donor")
dds <- DESeq(dds)
#Check this 
resultsNames(dds)
res <- results(dds, name = "origin_recipient_vs_donor", alpha = 0.05)
# pdf("MA_plot_before.pdf", width = 6, height = 6)
# plotMA(res)
# dev.off()
res <- lfcShrink(dds, coef = "origin_recipient_vs_donor", res = res, type = "apeglm")
# pdf("MA_plot_after.pdf", width = 6, height = 6)
# plotMA(res)
# dev.off()
res_tbl <- res %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  as_tibble() %>%
  filter(!is.na(padj)) %>%
  arrange(padj)
write.table(res_tbl,
            file = paste0("DESeq2_results_", gsub(" ", "_", ct), ".csv"),
            sep = ",", quote = FALSE, col.names = TRUE, row.names = FALSE)

# Add significance classification
res_tbl <- res_tbl %>%
  mutate(
    significance = case_when(
      padj < 0.05 & log2FoldChange > 1  ~ "Up in recipient",
      padj < 0.05 & log2FoldChange < -1 ~ "Up in donor",
      TRUE                              ~ "Not significant"
    )
  )


# Select top 20 genes for labeling
label_genes <- res_tbl %>%
  filter(abs(log2FoldChange) > 1) %>%
  arrange(padj) %>%
  slice_head(n = 20) %>%
  mutate(significance = case_when(
    padj < 0.05 & log2FoldChange > 1  ~ "Up in recipient",
    padj < 0.05 & log2FoldChange < -1 ~ "Up in donor",
    TRUE                              ~ "Not significant"
  ))

# Create volcano plot
volcano <- ggplot(res_tbl, aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
  geom_point(alpha = 0.7) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  ggrepel::geom_text_repel(
    data = label_genes,
    aes(label = gene),
    size = 3, max.overlaps = 15,
    show.legend = FALSE
  ) +
  scale_color_manual(values = c(
    "Up in recipient" = "#8B6A4F",  # Darker brown
    "Up in donor" = "#331F2C",
    "Not significant" = "gray"
  )) +
  labs(
    x = "Log2 Fold Change (recipient vs donor)",
    y = "-log10(padj)",
    color = "Direction"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.ticks = element_line(color = "black"),
    axis.line = element_line(color = "black")
  )


volcano

ggsave(paste0("Volcano_plot_", gsub(" ", "_", ct), ".pdf"), volcano, width = 6, height = 4)

# Select top 10 genes up in recipient
top_10_recipient <- res_tbl %>%
  filter(log2FoldChange > 1, padj < 0.05) %>%
  arrange(padj) %>%
  slice_head(n = 10) %>%
  pull(gene)

 
# Get normalized counts
  norm_counts <- counts(dds, normalized = TRUE) %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    filter(gene %in% top_10_recipient)
  
# Reshape for boxplot
  long_df <- norm_counts %>%
  pivot_longer(-gene, names_to = "Sample_name", values_to = "value") %>%
  mutate(Sample_name = gsub("\\.", "-", Sample_name),
           value = ifelse(value == 0, NA, value))
  
  # Get and join metadata
  metadata <- as.data.frame(colData(dds)) %>%
    rownames_to_column("Sample_name")
  
  plot_df <- left_join(long_df, metadata, by = "Sample_name") %>%
    filter(origin %in% c("donor", "recipient"))
  
# Boxplot
  top_10_plot <- ggplot(plot_df, aes(x = origin, y = value, color = origin)) +
    geom_boxplot(outlier.shape = NA, width = 0.6) +
    geom_jitter(width = 0.2, size = 2, alpha = 0.8) +
    scale_y_continuous(trans = "log10") +
    facet_wrap(~gene, scales = "free_y") +
    labs(title = paste("Top DE Genes in Recipient -", ct),
         x = "Origin", y = "Log10 normalized counts") +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 0))

  
ggsave(paste0("top10_recipient_DE_genes", gsub(" ", "_", ct), ".pdf"),
         plot = top_10_plot, width = 14, height = 10)
  
# Heatmap matrix
  mat <- norm_counts %>%
    column_to_rownames("gene") %>%
    as.matrix()
  mat <- mat[top_10_recipient, ]
  mat_scaled <- t(scale(t(mat)))  # z-score by gene
  
  # Plot heatmap
  pdf(paste0("heatmap_top10_", gsub(" ", "_", ct), ".pdf"), width = 8, height = 6)
  pheatmap::pheatmap(
    mat_scaled,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    annotation_col = as.data.frame(colData(dds)[, "origin", drop = FALSE]),
    show_rownames = TRUE,
    show_colnames = TRUE,
    main = paste("Top 10 DE Genes in", ct)
  )
  dev.off()
  


# ##### Loop version for all celltypes #####
# celltypes <- sort(unique(ex@meta.data$celltype))
# pb_list <- list()
# for (ct in celltypes) {
#   
#   # Subset cells to one celltype
#   seurat_ct <- subset(ex, subset=(celltype == ct))
#   
#   # Aggregate to get pseudobulk
#   bulk_ct <- AggregateExpression(
#     seurat_ct,
#     return.seurat = T,
#     assays = "RNA",
#     group.by = c("sample_id", "origin")
#   )
#   
#   # Fix sample_id formatting if needed
#   bulk_ct@meta.data$sample_id <- str_replace_all(bulk_ct@meta.data$sample_id, "_", "-")
#   
#   # Add number of cells per sample
#   n_cells <- seurat_ct@meta.data %>% 
#     dplyr::count(sample_id, origin) %>% 
#     dplyr::rename("n_cells"="n") %>%
#     mutate(orig.ident = paste0(str_replace_all(sample_id, "_", "-"), "_", origin)) %>%
#     dplyr::select(-c(origin, sample_id))
#   
#   meta_bulk_ct <- left_join(bulk_ct@meta.data, n_cells,  by ="orig.ident")
#   
#   head(meta_bulk_ct)
#   rownames(meta_bulk_ct) <- meta_bulk_ct$orig.ident
#   bulk_ct@meta.data <- meta_bulk_ct
#   pb_list[[ct]] <- bulk_ct
#   
# 
#   # Convert origin to factor for DESeq2
#   bulk_ct@meta.data$origin <- factor(bulk_ct@meta.data$origin)
#   
#   # Skip if not enough groups
#   if (length(unique(bulk_ct@meta.data$origin)) < 2) {
#     warning("Skipping ", ct, ": not enough groups")
#     next
#   }
#   
#   # Barplot
#   barplot <- ggplot(bulk_ct@meta.data, aes(x = sample_id, y = n_cells, fill = origin)) +
#     geom_bar(stat = "identity", color = "black") +
#     theme_classic() +
#     theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#     labs(title = paste("Cell Counts:", ct), x = "Sample", y = "Number of Cells") +
#     geom_text(aes(label = n_cells), vjust = -0.5)
#   ggsave(paste0("cell_counts_", gsub(" ", "_", ct), ".pdf"), barplot, width = 8, height = 5)
#   
#   # Fetch count matrix and run DESeq2
#   cluster_counts <- FetchData(bulk_ct, layer = "counts", vars = rownames(bulk_ct))
#   
#   # Check the order of samples is the same
#   all(rownames(cluster_counts)==rownames(bulk_ct@meta.data))
#   
#   # Create DESeq2 object
#   dds <- DESeqDataSetFromMatrix(t(cluster_counts), colData = bulk_ct@meta.data, design = ~ origin)
#   dds$origin <- relevel(dds$origin, ref="recipient")
#   dds <- DESeq(dds)
#   res <- results(dds, name = "origin_donor_vs_recipient", alpha = 0.05)
#   res <- lfcShrink(dds, coef = "origin_donor_vs_recipient", res = res, type = "apeglm")
#   
#   # Create DE result table
#   res_tbl <- res %>%
#     as.data.frame() %>%
#     rownames_to_column("gene") %>%
#     as_tibble() %>%
#     filter(!is.na(padj)) %>%
#     arrange(padj)
#   
#   # Save full results table
#   write.table(res_tbl,
#               file = paste0("DESeq2_results_", gsub(" ", "_", ct), ".csv"),
#               sep = ",", quote = FALSE, col.names = TRUE, row.names = FALSE)
#   
#   # Volcano plot
#   label_genes <- res_tbl %>%
#     filter(!is.na(padj)) %>%
#     filter(abs(log2FoldChange) > 1) %>%
#     arrange(padj) %>%
#     slice_head(n = 20)  # top 20 most significant (up/down)
#   
#   # Volcano plot with gene labels
#   volcano <- ggplot(res_tbl, aes(x = log2FoldChange, y = -log10(padj))) +
#     geom_point(alpha = 0.5) +
#     geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
#     geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
#     ggrepel::geom_text_repel(
#       data = label_genes,
#       aes(label = gene),
#       size = 3, max.overlaps = 15
#     ) +
#     labs(title = paste("Volcano plot -", ct),
#          x = "Log2 Fold Change",
#          y = "-log10(padj)") +
#     theme_minimal()
#   
#   ggsave(paste0("Volcano_plot_", gsub(" ", "_", ct), ".pdf"), barplot, width = 8, height = 5)
#   
#   # Top 10 DE genes up in recipient
#   top_10_recipient <- res_tbl %>%
#     filter(log2FoldChange > 1) %>%
#     slice_head(n = 10) %>%
#     pull(gene)
#   
#   if (length(top_10_recipient) > 1) {
#     norm_counts <- counts(dds, normalized = TRUE) %>%
#       as.data.frame() %>%
#       rownames_to_column("gene") %>%
#       filter(gene %in% top_10_recipient)
#     
#     # Reshape for plotting
#     long_df <- norm_counts %>%
#       pivot_longer(-gene, names_to = "Sample_name", values_to = "value") %>%
#       mutate(Sample_name = gsub("\\.", "-", Sample_name),
#              value = ifelse(value == 0, NA, value))  # handle log10(0)
#     
#     metadata <- as.data.frame(colData(dds)) %>%
#       rownames_to_column("Sample_name")
#     
#     plot_df <- left_join(long_df, metadata, by = "Sample_name") %>%
#       filter(origin %in% c("donor", "recipient"))
#     
#     # Boxplot
#     top_10_plot <- ggplot(plot_df, aes(x = origin, y = value, color = origin)) +
#       geom_boxplot(outlier.shape = NA, width = 0.6) +
#       geom_jitter(width = 0.2, size = 2, alpha = 0.8) +
#       scale_y_continuous(trans = "log10") +
#       facet_wrap(~gene, scales = "free_y") +
#       labs(title = paste("Top DE Genes in Recipient -", ct),
#            x = "Origin", y = "Log10 normalized counts") +
#       theme_minimal(base_size = 12) +
#       theme(plot.title = element_text(hjust = 0.5),
#             axis.text.x = element_text(angle = 0))
#     
#     ggsave(paste0("top10_recipient_DE_genes_", gsub(" ", "_", ct), ".pdf"),
#            plot = top_10_plot, width = 14, height = 10)
#     
#     # Heatmap
#     mat <- norm_counts %>%
#       column_to_rownames("gene") %>%
#       as.matrix()
#     mat <- mat[top_10_recipient, ]
#     mat_scaled <- t(scale(t(mat)))  # z-score by gene
#     
#     pdf(paste0("heatmap_top10_", gsub(" ", "_", ct), ".pdf"), width = 8, height = 6)
#     pheatmap(mat_scaled,
#              cluster_rows = TRUE,
#              cluster_cols = TRUE,
#              annotation_col = as.data.frame(colData(dds)[, "origin", drop = FALSE]),
#              show_rownames = TRUE, show_colnames = TRUE,
#              main = paste("Top 10 DE Genes in", ct))
#     dev.off()
#   } else {
#     message("Not enough DE genes to plot for ", ct)
#   }
# }

##end of the loop





