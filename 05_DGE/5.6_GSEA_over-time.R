# Nurefsan Sariipek and Peter van Galen, 250710

# Load libraries
library(tidyverse)
library(Seurat)
library(fgsea)

# Set working directory
repo_root <- system("git rev-parse --show-toplevel", intern = TRUE)
setwd(paste0(repo_root, "/05_DGE"))

# Clear environment variables
rm(list = ls())

# Load pathways for GSEA
pathways <- gmtPathways("c2.all.v2024.1.Hs.symbols.gmt") # used for paper
pathways <- gmtPathways("h.all.v2024.1.Hs.symbols.gmt") # out of curiosity

# Rank the DEG results from 5.5_DGE_tumorcells.R
de_results <- read_tsv("5.5_DGE_Pre-transplant_vs_Relapse.txt")
ranks <- as.numeric(de_results$log2FoldChange)
names(ranks) <- de_results$gene

# Run GSEA and order by p-value
gseaRes = fgsea(
  pathways = pathways,
  stats = ranks,
  minSize = 5,
  eps = 0.0,
  maxSize = 500,
  nPermSimple = 10000
)
gseaRes <- gseaRes %>% arrange(padj)

# Save results as tsv
data.table::fwrite(
  gseaRes,
  file = "5.6_GSEA_Pre-transplant_vs_Relapse.txt",
  sep = "\t",
  sep2 = c("", " ", "")
)

# Wrangle for plotting
df <- gseaRes %>%
  arrange(padj) %>%
  slice_head(n = 15) %>%
  filter(padj < 0.05) %>%
  dplyr::select(ID = pathway, padj, NES) %>%
  mutate(padj = -log10(padj)) %>%
  arrange(desc(NES)) %>%
  mutate(ID = factor(ID, levels = ID))

p1 <- ggplot(df, aes(x = NES, y = ID, fill = padj)) +
  geom_bar(stat = "identity") +
  geom_vline(xintercept = 0, color = "black") +
  geom_text(
    aes(label = ID),
    color = "black",
    hjust = ifelse(df$NES > 0, 1.05, -0.05)
  ) + # labels inside the bars
  scale_fill_gradient(
    low = "powderblue",
    high = "steelblue",
    name = "-log10(padj)"
  ) +
  coord_cartesian(clip = "off") +
  theme_minimal() +
  theme(
    aspect.ratio = 1,
    panel.grid = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(color = "black"),
    axis.ticks.x = element_line(color = "black"),
    axis.line.x = element_line(color = "black"),
    legend.position = "right",
    plot.margin = margin(t = 30, r = 20, b = 10, l = 10)
  ) +
  labs(x = "NES", y = NULL) +
  annotate(
    "text",
    x = -1.5,
    y = length(df$ID) + 1,
    label = "Up in pre-transplant",
    hjust = 0.5
  ) +
  annotate(
    "text",
    x = 1.5,
    y = length(df$ID) + 1,
    label = "Up in relapse",
    hjust = 0.5
  )

# View
p1

ggsave("5.6_GSEA_Relapse_pathways.pdf", height = 5, width = 5)
