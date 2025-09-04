# Nurefsan Sariipek and Peter van Galen, 250710
# Run the DGE analysis on recipient and donor HSPCs

# Load the needed libraries
library(tidyverse)
library(Seurat)
library(DESeq2)
library(apeglm)
library(pheatmap)

# Empty environment
rm(list = ls())

# Set working directory (Nurefsan)
setwd("~/hsct_immune_reconstitution/05_DGE/")
# For Peter:
setwd(
  "~/DropboxMGB/Projects/ImmuneEscapeTP53/hsct_immune_reconstitution/05_DGE/"
)

# Load the saved Seurat object
seu <- readRDS("../AuxiliaryFiles/250528_Seurat_complete.rds")

# Subset cells for DGE analysis
seu_subset <- seu %>%
  subset(
    sample_status == "remission" &
      celltype %in%
        c("HSC MPP", "MEP", "LMPP", "Cycling Progenitor", "Early GMP") &
      souporcell_origin %in% c("donor", "recipient")
  )

# Check cell numbers
seu_subset@meta.data %>%
  group_by(sample_id, cohort, souporcell_origin, timepoint) %>%
  dplyr::count() %>%
  pivot_wider(names_from = souporcell_origin, values_from = n) %>%
  print(n = 30)

# Subset further
seu_subset <- seu_subset %>%
  subset(cohort == "relapse" & timepoint %in% c(3, 5, 6))

# Also exclude P23 who does not have recipient HSPCs
seu_subset <- subset(seu_subset, sample_id != "P23_Rem1")
seu_subset$patient_id %>% unique %>% sort
# --> "For this analysis, we excluded the long-term remission cohort which had <10 persistent recipient HSPCs and P30-P33 who did not have remission samples at the 3-month time point. Thus, the following analysis of recipient vs. donor HSPCs during remission is based on six patients who developed relapse within 18 months."

# Check final cell numbers
n_cells <- seu_subset@meta.data %>%
  dplyr::count(patient_id, sample_id, cohort, souporcell_origin)
n_cells %>%
  pivot_wider(names_from = souporcell_origin, values_from = n) %>%
  adorn_totals(where = "row")
# --> "Genes that were upregulated in recipient HSPCs (n=533 cells) compared to their donor counterparts (n=540) included..."

# Visualize cell counts
n_cells %>%
  ggplot(aes(
    x = paste(sample_id, souporcell_origin),
    y = n,
    fill = souporcell_origin
  )) +
  geom_bar(stat = "identity", color = "black") +
  geom_text(aes(label = n), vjust = -0.5) +
  scale_fill_manual(values = c(recipient = "#E4C9B0", donor = "#4B3140")) +
  ylab("Number of cells") +
  theme_bw() +
  theme(
    aspect.ratio = 1,
    panel.grid = element_blank(),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_blank()
  )
ggsave("3.3.1_Cell_counts.pdf", width = 5, height = 5)

# Aggregate to get pseudobulk expression
seu_pseudobulk <- AggregateExpression(
  seu_subset,
  return.seurat = TRUE,
  assays = "RNA",
  group.by = c("sample_id", "souporcell_origin")
)

# Add metadata (number of cells per sample)
rownames(n_cells) <- paste0(
  gsub("_", "-", n_cells$sample_id),
  "_",
  n_cells$souporcell_origin
)
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
  design = ~ patient_id + souporcell_origin
)

# Specify reference condition
dds$souporcell_origin <- relevel(dds$souporcell_origin, ref = "donor")
dds <- DESeq(dds)

# Get results for the comparison of interest
resultsNames(dds)
res <- results(dds, name = "souporcell_origin_recipient_vs_donor", alpha = 0.05)

# Visualize
pdf("3.3.2_MA_plot_before_shrinking.pdf", width = 8, height = 6)
plotMA(res)
dev.off()

# Shrink results
res <- lfcShrink(
  dds,
  coef = "souporcell_origin_recipient_vs_donor",
  res = res,
  type = "apeglm"
)

pdf("3.3.3_MA_plot_after_shrinking.pdf", width = 6, height = 6)
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
  file = "3.3_DGE_Recipient_vs_Donor_HSPCs.txt",
)

# Add significance classification
res_tbl <- res_tbl %>%
  mutate(
    significance = case_when(
      padj < 0.05 & log2FoldChange > 1 ~ "Up in recipient",
      padj < 0.05 & log2FoldChange < -1 ~ "Up in donor",
      TRUE ~ "Not significant"
    )
  )

# Select top 20 genes for labeling
label_genes <- res_tbl %>%
  filter(abs(log2FoldChange) > 1) %>%
  arrange(padj) %>%
  slice_head(n = 20)

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
      "Up in recipient" = "#8B6A4F", #"#E4C9B0"
      "Up in donor" = "#331F2C", #"#4B3140"
      "Not significant" = "gray"
    )
  ) +
  labs(
    x = "Log2 fold change (recipient vs donor)",
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

ggsave("3.3.4_Volcano.pdf", volcano, width = 6, height = 4)

# Select top 10 genes up in recipient
top_10_recipient <- res_tbl %>%
  filter(log2FoldChange > 1, padj < 0.05) %>%
  arrange(padj) %>%
  slice_head(n = 10) %>%
  pull(gene)

# Get normalized counts
norm_counts <- counts(dds, normalized = T) %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  filter(gene %in% top_10_recipient)

# Reshape for boxplot
long_df <- norm_counts %>%
  pivot_longer(-gene, names_to = "Sample_name", values_to = "value")

# Get and join metadata
metadata <- as.data.frame(colData(dds)) %>%
  rownames_to_column("Sample_name")

plot_df <- left_join(long_df, metadata, by = "Sample_name") %>%
  filter(souporcell_origin %in% c("donor", "recipient"))

# Boxplot of top 10 upregulated genes in recipient HSPCs
top_10_boxplot <- plot_df %>%
  ggplot() +
  geom_boxplot(
    aes(x = souporcell_origin, y = value),
    outlier.shape = NA,
    width = 0.6
  ) +
  geom_jitter(
    aes(x = souporcell_origin, y = value, color = sample_id),
    width = 0.2,
    size = 2,
    alpha = 0.8
  ) +
  scale_y_continuous(trans = "log10") +
  facet_wrap(~gene, scales = "free_y", nrow = 2) +
  labs(
    title = paste("Top 10 genes in recipient vs. donor HSPCs"),
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
  paste0("3.3.5_Top10_recipient_boxplots.pdf"),
  plot = top_10_boxplot,
  width = 10,
  height = 4
)

# Heatmap matrix
mat <- norm_counts %>%
  column_to_rownames("gene") %>%
  as.matrix()
mat <- mat[top_10_recipient, ]
mat_scaled <- t(scale(t(mat))) # z-score by gene

# Plot heatmap
pdf("3.3.6_Top10_recipient_heatmap.pdf", width = 8, height = 6)
pheatmap::pheatmap(
  mat_scaled,
  cluster_rows = T,
  cluster_cols = T,
  annotation_col = as.data.frame(colData(dds)[,
    "souporcell_origin",
    drop = FALSE
  ]),
  annotation_colors = list(
    souporcell_origin = c("recipient" = "#E4C9B0", "donor" = "#4B3140")
  ),
  show_rownames = T,
  show_colnames = T,
  cellwidth = 20,
  cellheight = 20,
  main = paste("Top 10 genes in recipient vs. donor HSPCs")
)
dev.off()
