# Visualize GSEA systematically

library(dplyr)
library(ggpubr)
library(ComplexHeatmap)
library(circlize)


View(celltype_names)
t = read.table("~/Dropbox (Partners HealthCare)/ImmuneEscapeTP53/AnalysisNurefsan/DGE/posthsct/all_celltypes.total_GSEA_result.tsv", sep = "\t", header = T)
View(t)

celltypes = c("Monocytes", "CD4 Memory", "Mid Erythroids", "CD8 Naïve", "CD4 Naïve", "CD8 Effector Memory", "NK T cells", "Late Erythroids", "Early Erythroids","CD8 Central Memory") 
t$celltype[t$celltype == "CD8 Effector Memory"] = "CD8 Effector"
t$celltype[t$celltype == "CD8 Central Memory"] = "CD8 Memory"

t$celltype = factor(t$celltype, levels = celltypes)
sum(is.na(t$celltype))



t %>% 
  group_by(ID) %>%
  summarize(nhits = sum(p.adjust<0.05), median_NES = median(NES)) %>%
  arrange(desc(nhits)) %>%
  filter( grepl("HALLMARK|KEGG|GOBP", ID)) %>%
  ggplot(aes(x=nhits)) + geom_histogram(binwidth=1) + scale_y_log10() + theme_pubr() + ggtitle("Immunologic gene sets")

relevant_ids = t %>% 
  group_by(ID) %>%
  summarize(nhits = sum(p.adjust<0.05), median_NES = median(NES)) %>%
  arrange(desc(nhits)) %>%
  filter( grepl("HALLMARK|KEGG|GOBP", ID)) %>%
  filter(nhits > 0) %>%
  pull(ID) %>%
  unique()


p = t %>%
  filter(ID %in% grep("GOBP", relevant_ids, value = T)) %>%
  dplyr::select(ID, NES, p.adjust, celltype) %>%
  mutate(
    label = case_when(
      p.adjust > 0.1 ~ "",
      p.adjust > 0.01 ~ "*",
      p.adjust > 0.001 ~ "**",
      !is.na(p.adjust) ~ "***",
      TRUE ~ NA_character_
    )) %>%
  mutate(group = case_when(
    celltype %in% c("Early Erythroids", "Late Erythroids", "Mid Erythroids") ~ "Erythroid",
    celltype == "Monocytes" ~ "Myeloid",
    celltype %in% c("CD4 Memory", "CD8 Naïve", "CD4 Naïve", "CD8 Effector Memory","CD8 Central Memory",  "NK T cells") ~ "T cells")) %>%
  ggplot(aes(y=NES, x = celltype, fill = group)) + 
  geom_bar(stat = "identity", position = "dodge", width = 0.75, linewidth = 1.5) + 
  facet_wrap(ID~., ncol=3) +
  geom_text(aes(label = label, y = ifelse(NES>0, NES+0.1, NES-0.3))) + 
  theme_pubr(base_size = 10) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.99)) + xlab("") +
  scale_fill_manual(values = c("#dddddd", "#d68e7e", "#f0c15d", "#9fd1cf",  "#d6bce3"), breaks = c     ( "Erythroid", "Myeloid", "T cells"), name = "Lineage") + 
  scale_color_manual(values = c("grey30", "brown3"))

p

save_plot("top_hallmark_pathways_barplot.pdf", p, base_height = 20, base_width = 30, limitsize = FALSE) 



