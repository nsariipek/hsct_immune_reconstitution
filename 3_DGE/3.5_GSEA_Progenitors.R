# Load the libraries
library(fgsea)
library(dplyr)
library(Seurat)

# Empty environment
rm(list=ls())

# Set working directory
setwd("~/TP53_ImmuneEscape/3_DGE/")

# Define the pathways, which can be found in the google bucket
c2_pathways = gmtPathways("~/c2.all.v2024.1.Hs.symbols.gmt")
hallmark_pathways = gmtPathways("~/h.all.v2024.1.Hs.symbols.gmt")

# Define the DEG results from 3.4 that is saved in the folder
de_results = read.table("DESeq2_results_Merged_Progenitors.csv", header = T, sep = ",")
ranks = as.numeric(de_results$log2FoldChange)
names(ranks) = de_results$gene


gseaRes = fgsea(pathways = c2_pathways, 
                stats    = ranks,
                minSize  = 5,
                eps      = 0.0,
                maxSize  = 500,
                nPermSimple = 10000)

gseaRes = gseaRes %>% arrange(padj)

data.table::fwrite(gseaRes, file="DESeq2_GSEA_hallmarkpathways_results.txt", sep="\t", sep2=c("", " ", ""))

df = gseaRes %>%
  filter(padj<0.05) %>%
  arrange(padj) %>%
  slice_head(n=25) %>%
  dplyr::select(ID=pathway, padj, NES) %>%
  mutate(padj = -log10(padj)) %>%
  arrange((NES)) %>%
  mutate(ID = factor(ID, levels = ID))

p1 = df %>%
  ggplot(aes(x=NES, y=ID, fill = padj)) +
  geom_bar(stat = "identity") +
  theme_pubr(base_size = 8) +
  scale_fill_gradient(low = "#D2B48C", high = "#8B4513", trans="log", breaks = scales::pretty_breaks(), name =  expression(-log[10] ~ P[adj])) +
  theme(axis.ticks = element_line(color = "black"),
        axis.line = element_line(linewidth=0.3), legend.direction = "vertical",
        legend.position = "right",  # Moves the legend above the plot
        legend.justification = "left",
        legend.key.width = unit(0.4, "cm"),
        legend.key.height = unit(0.3, "cm"))
p1

ggsave("gsea_plot_c2_sig.pdf", height = 6, width = 8)
dev.off()

