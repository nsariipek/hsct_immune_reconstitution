# Peter van Galen, 251129
# Generate informative sina- and barplots

# First question by Andy (Nov '25): "Did you get any hits for CathepsinG in the transplant MRD study? We should also look in our ReDeeM AML data. There is a trial of this T cell engager already happening in R/ RAML, and the BTC group is look at AML MRD." - https://aacrjournals.org/cancerres/article/84/6_Supplement/1236/739282/Abstract-1236-Characterization-of-CBX-250-a-first

# Load libraries
library(tidyverse)
library(Seurat)
library(ggforce)
#library(cowplot)

# Set working directory
repo_root <- system("git rev-parse --show-toplevel", intern = TRUE)
setwd(paste0(repo_root, "/AnalysisPeter"))

# Clear environment variables
rm(list = ls())

# Load and normalize data
seu <- readRDS("../AuxiliaryFiles/250528_Seurat_complete.rds")
seu <- NormalizeData(seu)

# Plot single gene ------------------------------------------------------------

# Generate a metadata tibble and add gene expression (choose one gene)
metadata <- as_tibble(seu@meta.data, rownames = "cell")
mygene <- "CTSG"
mygene <- "PRAME"
metadata$mygene <- LayerData(seu, layer = "data")[mygene, ]


# Expression in donor/recipient HSPCs -----------------------------------------

# Get P-value from original analysis (pseudobulk for "HSC MPP", "MEP", "LMPP", "Cycling Progenitor", "Early GMP")
res_tbl <- read_tsv("../05_DGE/5.3_DGE_Recipient_vs_Donor_HSPCs.txt")
res_tbl %>% filter(gene == mygene)

# Visualize
metadata %>%
  filter(
    sample_status == "remission" &
      celltype %in%
        c("HSC MPP", "MEP", "LMPP", "Cycling Progenitor", "Early GMP") &
      souporcell_origin %in% c("donor", "recipient")
  ) %>%
  # Reorder facets
  mutate(
    souporcell_origin = factor(
      souporcell_origin,
      levels = c("donor", "recipient")
    )
  ) %>%
  ggplot(aes(x = celltype, y = mygene, color = sample_id)) +
  geom_sina(aes(group = celltype), size = 0.1, scale = "width") +
  geom_violin(aes(group = celltype), fill = NA, scale = "width") +
  facet_wrap(~souporcell_origin) +
  ggtitle(paste0(
    mygene,
    " in recipient vs. donor HSPCs (padj = ",
    format(filter(res_tbl, gene == mygene)$padj, digits = 3),
    ", DESeq2 pseudobulked)"
  )) +
  ylab(paste(mygene, "expression")) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )
#ggsave(paste0("251129_", mygene, "_expression_HSPCs.pdf"))

# Expression in donor/recipient myeloid cells ---------------------------------

# fmt: skip
cell_subset <- c("HSC MPP", "MEP", "EoBasoMast Precursor", "Megakaryocyte Precursor", "Early Erythroid",
"Late Erythroid", "LMPP", "Cycling Progenitor" ,"Early GMP", "Late GMP", "Pro-Monocyte",
"CD14 Mono", "CD16 Mono", "cDC")

# Calculate p-values for each cell type (recipient vs. donor) using Wilcox test. Note! This is very permissive
pvalue_results <- metadata %>%
  filter(
    sample_status == "remission" &
      celltype %in% cell_subset &
      souporcell_origin %in% c("donor", "recipient")
  ) %>%
  group_by(celltype) %>%
  summarize(
    n_donor = sum(souporcell_origin == "donor"),
    n_recipient = sum(souporcell_origin == "recipient"),
    mean_donor = mean(mygene[souporcell_origin == "donor"]),
    mean_recipient = mean(mygene[souporcell_origin == "recipient"]),
    pvalue = if (n_donor > 0 & n_recipient > 0) {
      wilcox.test(
        mygene[souporcell_origin == "donor"],
        mygene[souporcell_origin == "recipient"]
      )$p.value
    } else {
      NA_real_
    },
    .groups = "drop"
  ) %>%
  mutate(padj = p.adjust(pvalue, method = "BH")) %>%
  mutate(
    sig_label = case_when(
      padj < 0.001 ~ "***",
      padj < 0.01 ~ "**",
      padj < 0.05 ~ "*",
      TRUE ~ ""
    )
  )

# View results
pvalue_results

# Expression in donor/recipient myeloid cells
metadata %>%
  filter(
    sample_status == "remission" &
      celltype %in% cell_subset,
    souporcell_origin %in% c("donor", "recipient")
  ) %>%
  ggplot(aes(x = celltype, y = mygene, color = sample_id)) +
  geom_sina(aes(group = celltype), size = 0.1, scale = "width") +
  geom_violin(aes(group = celltype), fill = NA, scale = "width") +
  geom_text(
    data = pvalue_results %>% mutate(souporcell_origin = "recipient"),
    aes(x = celltype, y = Inf, label = sig_label),
    vjust = 1.5,
    hjust = 0.5,
    color = "black",
    inherit.aes = FALSE
  ) +
  geom_label(
    data = pvalue_results %>%
      pivot_longer(
        cols = c(mean_donor, mean_recipient),
        names_to = "souporcell_origin",
        values_to = "mean"
      ) %>%
      mutate(
        souporcell_origin = ifelse(
          souporcell_origin == "mean_donor",
          "donor",
          "recipient"
        )
      ),
    aes(x = celltype, y = mean, label = round(mean, 2)),
    vjust = -0.5,
    hjust = 0.5,
    color = "black",
    inherit.aes = FALSE
  ) +
  facet_wrap(~souporcell_origin) +
  ggtitle(paste0(
    "Mean ",
    mygene,
    " in recipient vs. donor myeloid cells (P-values from Wilcoxon test, very lenient)"
  )) +
  ylab(paste(mygene, "expression")) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )
#ggsave(paste0("251129_", mygene, "_expression_myeloid.pdf"), width = 10)

# Heatmap of all cell types, separate by patient ------------------------------

# This requires >378 Gb RAM and 20 minutes (such as Google VM n2-highmem-64)
seu <- ScaleData(seu)

# Generate a metadata tibble and add gene expression (choose one gene)
metadata_h <- as_tibble(seu@meta.data, rownames = "cell")
mygene_h <- "PRAME"
metadata_h$mygene_h <- LayerData(seu, layer = "scale.data")[mygene_h, ]


metadata_h %>%
  filter(
    sample_status == "remission",
    celltype %in% levels(celltype),
    souporcell_origin %in% c("donor", "recipient")
  ) %>%
  mutate(
    souporcell_origin = factor(
      souporcell_origin,
      levels = c("donor", "recipient")
    )
  ) %>%
  group_by(sample_id, celltype, souporcell_origin) %>%
  summarize(n = n(), mean_mygene_h = mean(mygene_h), .groups = "drop") %>%
  # fmt: skip
  complete(sample_id, celltype, souporcell_origin,
    fill = list(n = 0, mean_mygene_h = NA)) %>%
  # Remove samples where all celltypes have 0 cells
  group_by(sample_id, souporcell_origin) %>%
  filter(sum(n) > 0) %>%
  ungroup() %>%
  mutate(outline = if_else(n < 10 & n > 0, "grey", NA_character_)) %>%
  ggplot(aes(x = celltype, y = sample_id, fill = mean_mygene_h)) +
  geom_tile(aes(color = outline), linewidth = 0.5) +
  scale_color_identity() +
  facet_wrap(~souporcell_origin) +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    name = "Scaled\nExpression",
    na.value = "grey80"
  ) +
  labs(
    title = paste0(
      "Mean scaled ",
      mygene_h,
      " expression at remission (grey fill indicates 0 cells, grey box indicates <10 cells)"
    )
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )
ggsave(paste0("251129_", mygene_h, "_heatmap.pdf"), width = 15)
