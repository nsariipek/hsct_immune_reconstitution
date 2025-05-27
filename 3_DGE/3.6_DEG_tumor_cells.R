# Nurefsan Sariipek, 250522
# DGE comparison in between tumor cells
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

# Set working directory (Nurefsan)
setwd("~/TP53_ImmuneEscape/3_DGE/")

# Load the seurat object that has Numbat results from 8.4
seu_combined <- readRDS("~/250505_numbat_combined_seurat.rds")

# Set sample status order
seu_combined$sample_status <- factor(
  seu_combined$sample_status, levels = c("pre-transplant", "remission", "relapse"))

seu_combined$compartment_opt <- factor(seu_combined$compartment_opt,levels = c("normal", "tumor"))

# Normalize
seu_combined <- NormalizeData(seu_combined)
seu_combined <- FindVariableFeatures(seu_combined)
# Scale the data
seu_combined <- ScaleData(seu_combined)

# Select the group you want to run analysis
seu_subset <- seu_combined %>% subset(compartment_opt == "tumor" & sample_status %in% c("pre-transplant","relapse")
                                  

# & celltype %in% c("Late Erythroid", "Early Erythroid", "MEP","CD14 Mono", "CD16 Mono", "Pro-Monocyte", "Early GMP", "Late GMP","HSC MPP", "LMPP", "Cycling Progenitor", "Early GMP", "Late GMP", "MEP", "EoBasoMast Precursor", "Megakaryocyte Precursor","cDC", "pDC")
     #      & !sample_id %in% c("P23_Rem2","P31_Rem")
)
                                    
# Combine these cell types into one label
seu_subset$celltype_merged <- "Merged_cells"

# Check that how many cells are per pseudo-sample
seu_subset@meta.data %>%
  group_by(sample_id, celltype) %>%
  summarize(n=n()) %>% 
  pivot_wider(names_from = "celltype", values_from = "n") %>% View()

# Use all cells in seu_subset (which are now labeled Merged_Cells)
seurat_ct <- seu_subset
ct <- "Merged_cells"

# Aggregate to get pseudobulk
bulk_ct <- AggregateExpression(
  seurat_ct,
  return.seurat = TRUE,
  assays = "RNA",
  group.by = c("sample_id","sample_status")
)

# Format sample_id
bulk_ct@meta.data$sample_id <- str_replace_all(bulk_ct@meta.data$sample_id, "_", "-")

# Add number of cells per sample
n_cells <- seurat_ct@meta.data %>%
  dplyr::count(sample_id, sample_status) %>%
  dplyr::rename("n_cells" = "n") %>%
  mutate(orig.ident = paste0(str_replace_all(sample_id, "_", "-"), "_", sample_status)) %>%
  dplyr::select(-c(sample_status, sample_id))

meta_bulk_ct <- left_join(bulk_ct@meta.data, n_cells, by = "orig.ident")
rownames(meta_bulk_ct) <- meta_bulk_ct$orig.ident
bulk_ct@meta.data <- meta_bulk_ct

# Exclude unpaired samples
bulk_ct@meta.data$sample_status <- factor(bulk_ct@meta.data$sample_status)
if (length(unique(bulk_ct@meta.data$sample_status)) < 2) stop("Not enough groups to compare.")

# Check the cell counts per group
barplot <- ggplot(bulk_ct@meta.data, aes(x = sample_id, y = n_cells, fill = sample_status)) +
  geom_bar(stat = "identity", color = "black") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = paste("Cell Counts:", ct), x = "Sample", y = "Number of Cells") +
  geom_text(aes(label = n_cells), vjust = -0.5)
barplot
#ggsave(paste0("cell_counts_progenitors_", gsub(" ", "_", ct), ".pdf"), barplot, width = 8, height = 5)

#Run DESeq 
cluster_counts <- FetchData(bulk_ct, layer = "counts", vars = rownames(bulk_ct))

dds <- DESeqDataSetFromMatrix(t(cluster_counts), colData = bulk_ct@meta.data, design = ~sample_status) # we decided ~ sample_id + origin is better 
# This line put donor as a reference, but check line 136
dds$sample_status <- relevel(dds$sample_status, ref = "pre-transplant")
dds <- DESeq(dds)
#Check this 
resultsNames(dds)
res <- results(dds, name = "sample_status_relapse_vs_pre.transplant", alpha = 0.05)
# pdf("MA_plot_before.pdf", width = 6, height = 6)
#plotMA(res)
# dev.off()
res <- lfcShrink(dds, coef = "sample_status_relapse_vs_pre.transplant", res = res, type = "apeglm")
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
            file = paste0("DESeq2_results_sample_status_relapse_vs_pre.transplant", gsub(" ", "_", ct), ".csv"),
            sep = ",", quote = FALSE, col.names = TRUE, row.names = FALSE)

# Add significance classification
res_tbl <- res_tbl %>%
  mutate(
    significance = case_when(
      padj < 0.05 & log2FoldChange > 1  ~ "Up in relapse",
      padj < 0.05 & log2FoldChange < -1 ~ "Up in pre-transplant",
      TRUE                              ~ "Not significant"
    )
  )


# Select top 20 genes for labeling
label_genes <- res_tbl %>%
  filter(abs(log2FoldChange) > 1) %>%
  arrange(padj) %>%
  slice_head(n = 20) %>%
  mutate(
    significance = case_when(
      padj < 0.05 & log2FoldChange > 1  ~ "Up in relapse",
      padj < 0.05 & log2FoldChange < -1 ~ "Up in pre-transplant",
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
    "Up in pre-transplant" = "blue" ,  
    "Up in relapse" = "red",
    "Not significant" = "gray"
  )) +

  labs(
    x = "Log2 Fold Change (pre-transplant vs remission)",
    y = "-log10(padj)",
    color = "Direction"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.ticks = element_line(color = "black"),
    axis.line = element_line(color = "black")
  )


volcano

ggsave(paste0("Volcano_plot_sample_status_relapse_vs_pre.transplant_", gsub(" ", "_", ct), ".pdf"), volcano, width = 6, height = 4)

# Select top 10 genes up in relapse
top_10_relapse <- res_tbl %>%
  filter(log2FoldChange > 1, padj < 0.05) %>%
  arrange(padj) %>%
  slice_head(n = 10) %>%
  pull(gene)


# Get normalized counts
norm_counts <- counts(dds, normalized = TRUE) %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  filter(gene %in% top_10_relapse)

# Reshape for boxplot
long_df <- norm_counts %>%
  pivot_longer(-gene, names_to = "Sample_name", values_to = "value") %>%
  mutate(Sample_name = gsub("\\.", "-", Sample_name),
         value = ifelse(value == 0, NA, value))

# Get and join metadata
metadata <- as.data.frame(colData(dds)) %>%
  rownames_to_column("Sample_name")

plot_df <- left_join(long_df, metadata, by = "Sample_name") %>%
  filter(sample_status %in% c("pre-transplant", "relapse"))

# Boxplot
top_10_plot <- ggplot(plot_df, aes(x = sample_status, y = value, color = sample_status)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.8) +
  scale_y_continuous(trans = "log10") +
  facet_wrap(~gene, scales = "free_y") +
  labs(title = paste("Top DE Genes in Relapse-", ct),
       x = "sample_status", y = "Log10 normalized counts") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 0))
top_10_plot

ggsave(paste0("top10_relapse_DE_genes", gsub(" ", "_", ct), ".pdf"),
       plot = top_10_plot, width = 14, height = 10)

# Heatmap matrix
mat <- norm_counts %>%
  column_to_rownames("gene") %>%
  as.matrix()
mat <- mat[top_10_relapse, ]
mat_scaled <- t(scale(t(mat)))  # z-score by gene

# Plot heatmap
pdf(paste0("heatmap_top10_relapse_Vspre_transplant", gsub(" ", "_", ct), ".pdf"), width = 8, height = 6)
pheatmap::pheatmap(
  mat_scaled,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  annotation_col = as.data.frame(colData(dds)[, "sample_status", drop = FALSE]),
  show_rownames = TRUE,
  show_colnames = TRUE,
  main = paste("Top 10 DE Genes in", ct)
)
dev.off()
