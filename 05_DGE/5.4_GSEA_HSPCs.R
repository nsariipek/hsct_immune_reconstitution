# Nurefsan Sariipek and Peter van Galen, 250711

# Load the libraries
library(tidyverse)
library(Seurat)
library(fgsea)

# Empty environment
rm(list = ls())

# Set working directory
setwd("~/hsct_immune_reconstitution/05_DGE/")

# For Peter
setwd(
  "~/DropboxMGB/Projects/ImmuneEscapeTP53/hsct_immune_reconstitution/05_DGE"
)

# Load pathways for GSEA
hallmark_pathways <- gmtPathways("h.all.v2024.1.Hs.symbols.gmt")

# Rank the DGE results from 5.3_DGE_progenitors.R
de_results <- read_tsv("5.3_DGE_Recipient_vs_Donor_HSPCs.tsv")
ranks <- as.numeric(de_results$log2FoldChange)
names(ranks) <- de_results$gene

# Run GSEA and order by p-value
gseaRes <- fgsea(
  pathways = hallmark_pathways,
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
  file = "5.4_GSEA_Recipient_vs_Donor_HSPCs.txt",
  sep = "\t",
  sep2 = c("", " ", "")
)

# Wrangle for plotting
df <- gseaRes %>%
  arrange(padj) %>%
  filter(padj < 0.05) %>%
  dplyr::select(ID = pathway, padj, NES) %>%
  mutate(padj = -log10(padj)) %>%
  arrange(desc(NES)) %>%
  #  mutate(ID = gsub("HALLMARK_", "", ID)) %>%
  mutate(ID = factor(ID, levels = ID))

# Plot
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
    label = "Up in donor HSPCs",
    hjust = 0.5
  ) +
  annotate(
    "text",
    x = 1.5,
    y = length(df$ID) + 1,
    label = "Up in recipient HSPCs",
    hjust = 0.5
  )

# View
p1

ggsave("5.4_GSEA_HSPCs_Hallmark.pdf", height = 4, width = 8)
