# Load the libraries
library(fgsea)
library(dplyr)
library(Seurat)
library(colorspace)
library(ggnewscale)

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

gseaRes = fgsea(pathways = hallmark_pathways, 
                stats    = ranks,
                minSize  = 5,
                eps      = 0.0,
                maxSize  = 500,
                nPermSimple = 10000)

gseaRes = gseaRes %>% arrange(padj)

data.table::fwrite(gseaRes, file="DESeq2_GSEA_hallmarkpathways_results.txt", sep="\t", sep2=c("", " ", ""))


df = gseaRes %>%
  arrange(padj) %>%
  filter(padj<0.05) %>%
  slice_head(n=20) %>%
  dplyr::select(ID=pathway, padj, NES) %>%
  mutate(padj = -log10(padj)) %>%
  arrange(desc(NES)) %>%
  mutate(ID = factor(ID, levels = ID))

p1 <- ggplot(df, aes(x = NES, y = ID, fill = padj)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.3) +
  geom_text(aes(label = ID),
            color = "black", size = 2.3,
            hjust = ifelse(df$NES > 0, 1.05, -0.05)) +  # labels inside the bars
  scale_fill_gradient(
    low = "#E4C9B0",   # light beige
    high = "#4B3140",  # dark brown
    name = expression(-log[10] ~ P[adj])
  ) +
  scale_x_continuous(limits = c(-2.5, 2.5), expand = c(0, 0)) +
  theme_void(base_size = 8) +
  theme(
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 7),
    axis.ticks.x = element_line(color = "black", linewidth = 0.3),
    axis.line.x = element_line(color = "black", linewidth = 0.3),
    legend.position = "right",
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 6),
    legend.key.width = unit(1.0, "cm"),     
    legend.key.height = unit(0.3, "cm"),  
    plot.margin = margin(t = 30, r = 20, b = 10, l = 10) 
  )+
  labs(x = "NES", y = NULL) +
  annotate("text", x = -1, y = length(df$ID) + 0.8, label = "← UP in Donor Cells", size = 3, hjust = 0) +
  annotate("text", x =  1, y = length(df$ID) + 0.8, label = "UP in Recipient Cells →", size = 3, hjust = 1)




p1
ggsave("gsea_plot_c2_sig.pdf", height = 8, width = 8)


