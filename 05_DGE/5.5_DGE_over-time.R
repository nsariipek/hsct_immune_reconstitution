# Nurefsan Sariipek and Peter van Galen, 250706
# DGE comparison in between tumor cells

# Load the needed libraries
library(tidyverse)
library(Seurat)
library(DESeq2)
library(apeglm)
library(pheatmap)
library(janitor)

# Set working directory (Nurefsan)
setwd("~/hsct_immune_reconstitution/05_DGE/")

# For Peter
# fmt: skip
setwd("~/DropboxMGB/Projects/ImmuneEscapeTP53/hsct_immune_reconstitution/05_DGE")

# Remove environment variables
rm(list = ls())

# Load the saved Seurat object
seu <- readRDS("../AuxiliaryFiles/250528_Seurat_complete.rds")

#  Subset malignant cells for DGE analysis
seu_subset <- seu %>%
  subset(
    numbat_compartment == "tumor" &
      sample_status %in% c("pre-transplant", "remission", "relapse"),
  )

# Check cell numbers
n_cells <- seu_subset@meta.data %>%
  dplyr::count(patient_id, sample_status)

# Visualize cell numbers
n_cells %>%
  ggplot(aes(x = sample_status, y = n, fill = sample_status)) +
  geom_bar(stat = "identity", color = "black") +
  geom_text(aes(label = n), vjust = -0.5) +
  scale_fill_manual(
    values = c(
      "pre-transplant" = "#A3BFD9",
      "remission" = "#F6E06E",
      "relapse" = "#8B0000"
    )
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  facet_wrap(~patient_id) +
  theme_bw() +
  ylab("Number of cells") +
  theme(
    aspect.ratio = 1,
    panel.grid = element_blank(),
    axis.title.x = element_blank(),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
ggsave("5.5.1_Cell_counts.pdf", width = 7, height = 4)

# Exclude P22 and P33 who have only one timepoint, and exclude remission due to insufficient cell numbers
seu_subset <- subset(
  seu_subset,
  !patient_id %in% c("P22", "P33") & sample_status != "remission"
)

# Check final cell numbers
n_cells <- seu_subset@meta.data %>%
  dplyr::count(patient_id, sample_status)
n_cells %>%
  pivot_wider(names_from = sample_status, values_from = n) %>%
  adorn_totals(where = c("row"))
# --> In four patients, we detected sufficient malignant cells pre-transplant and at relapse (n=2,911 and n=31,763, Supplementary Figure 6A) to evaluate malignant cell evolution over time.

# Aggregate to get pseudobulk expression
seu_pseudobulk <- AggregateExpression(
  seu_subset,
  return.seurat = TRUE,
  assays = "RNA",
  group.by = c("patient_id", "sample_status")
)

# Add metadata (number of cells per sample)
rownames(n_cells) <- paste0(n_cells$patient_id, "_", n_cells$sample_status)
seu_pseudobulk <- AddMetaData(seu_pseudobulk, n_cells)
seu_pseudobulk@meta.data

# Run DESeq
cluster_counts <- FetchData(
  seu_pseudobulk,
  layer = "counts",
  vars = rownames(seu_pseudobulk)
)
dds <- DESeqDataSetFromMatrix(
  t(cluster_counts),
  colData = seu_pseudobulk@meta.data,
  design = ~ patient_id + sample_status
)

# Specify reference condition
dds$sample_status <- relevel(dds$sample_status, ref = "pre-transplant")
dds <- DESeq(dds)

# Get results
resultsNames(dds)
res <- results(
  dds,
  name = "sample_status_relapse_vs_pre.transplant",
  alpha = 0.05
)

# Visualize
pdf("5.5.2_MA_plot_before_shrinking.pdf", width = 8, height = 6)
plotMA(res)
dev.off()

# Shrink results
res <- lfcShrink(
  dds,
  coef = "sample_status_relapse_vs_pre.transplant",
  res = res,
  type = "apeglm"
)

pdf("5.5.3_MA_plot_after_shrinking.pdf", width = 6, height = 6)
plotMA(res)
dev.off()

# Wrangle and save results
res_tbl <- res %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  as_tibble() %>%
  filter(!is.na(padj)) %>%
  arrange(padj)
write_tsv(
  mutate(
    res_tbl,
    across(where(is.numeric), ~ format(.x, scientific = FALSE, digits = 15))
  ),
  file = "5.5_DGE_Pre-transplant_vs_Relapse.txt"
)

# Add significance classification
res_tbl <- res_tbl %>%
  mutate(
    significance = case_when(
      padj < 0.05 & log2FoldChange > 1 ~ "Up in relapse",
      padj < 0.05 & log2FoldChange < -1 ~ "Up in pre-transplant",
      TRUE ~ "Not significant"
    )
  )

# Select top genes for labeling
label_genes <- res_tbl %>%
  filter(abs(log2FoldChange) > 1) %>%
  arrange(padj) %>%
  slice_head(n = 20)
# Add some labels more on the right
label_genes <- rbind(
  label_genes,
  res_tbl %>% filter(log2FoldChange > 6),
  res_tbl %>% filter(log2FoldChange > 1 & padj < 10^-7)
) %>%
  unique()

# Create volcano plot
volcano <- res_tbl %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
  geom_point(alpha = 0.7) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  ggrepel::geom_text_repel(
    data = label_genes,
    aes(label = gene),
    size = 3,
    max.overlaps = 15,
    show.legend = FALSE
  ) +
  scale_color_manual(
    values = c(
      "Up in pre-transplant" = "#A3BFD9",
      "Up in relapse" = "#8B0000",
      "Not significant" = "gray"
    )
  ) +
  labs(
    x = "Log2 fold change (relapse vs. pre-transplant)",
    y = "-log10(padj)",
    color = "Direction"
  ) +
  theme_bw() +
  theme(
    aspect.ratio = 1,
    panel.grid = element_blank(),
    axis.line = element_blank(),
    axis.text = element_text(color = "black"),
    axis.ticks = element_line(color = "black"),
  )

# View
volcano

# Save
ggsave("5.5.4_Volcano.pdf", volcano, width = 6, height = 4)

# Select top 10 genes upregulated in relapse
top_10_relapse <- res_tbl %>%
  filter(log2FoldChange > 1, padj < 0.05) %>%
  arrange(padj) %>%
  slice_head(n = 10) %>%
  pull(gene)

# Get normalized counts
norm_counts <- counts(dds, normalized = T) %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  filter(gene %in% top_10_relapse)

# Reshape for boxplot
long_df <- norm_counts %>%
  pivot_longer(-gene, names_to = "Sample_name", values_to = "value")

# Get and join metadata
metadata <- as.data.frame(colData(dds)) %>%
  rownames_to_column("Sample_name")

plot_df <- left_join(long_df, metadata, by = "Sample_name") %>%
  filter(sample_status %in% c("pre-transplant", "relapse"))

# Boxplot of top 10 upregulated genes with relapse
top_10_boxplot <- plot_df %>%
  ggplot() +
  geom_boxplot(
    aes(x = sample_status, y = value),
    outlier.shape = NA,
    width = 0.6
  ) +
  geom_jitter(
    aes(x = sample_status, y = value, color = patient_id),
    width = 0.2,
    size = 2,
    alpha = 0.8
  ) +
  scale_y_continuous(trans = "log10") +
  facet_wrap(~gene, scales = "free_y", nrow = 2) +
  labs(
    title = paste("Top 10 genes in relapse vs. pre-transplant"),
    y = "Log10 normalized counts"
  ) +
  theme_bw() +
  theme(
    aspect.ratio = 1,
    panel.grid = element_blank(),
    axis.line = element_blank(),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.ticks = element_line(color = "black"),
    plot.title = element_text(hjust = 0.5)
  )

# View
top_10_boxplot

ggsave(
  paste0("5.5.5_Top10_relapse_boxplots.pdf"),
  plot = top_10_boxplot,
  width = 8,
  height = 6
)

# Heatmap matrix
mat <- norm_counts %>%
  column_to_rownames("gene") %>%
  as.matrix()
mat_scaled <- t(scale(t(mat))) # z-score by gene

# Plot heatmap
pdf("5.5.6_Top10_relapse_heatmap.pdf", width = 7, height = 6)
pheatmap::pheatmap(
  mat_scaled,
  cluster_rows = T,
  cluster_cols = T,
  annotation_col = as.data.frame(colData(dds)[, "sample_status", drop = FALSE]),
  annotation_colors = list(
    sample_status = c("pre-transplant" = "#A3BFD9", "relapse" = "#8B0000")
  ),
  show_rownames = T,
  show_colnames = T,
  cellwidth = 20,
  cellheight = 20,
  main = paste("Top 10 genes in relapse vs. pre-transplant")
)
dev.off()
