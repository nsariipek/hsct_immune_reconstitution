# Visualize GSEA systematically

library(dplyr)
library(ggpubr)
library(ComplexHeatmap)
library(circlize)
library(tidyr)

# Empty environment
rm(list=ls())

# For Nurefsan:
my_wd <- "/Users/dz855/Dropbox (Partners HealthCare)/ImmuneEscapeTP53/AnalysisNurefsan/"

t = read.table(paste0(my_wd,"DGE/only_Tcells-nonsig/GSEA_donorcells_result.tsv"), sep = "\t", header = T)
View(t)


celltypes = c("CD4 Naïve", "CD4 Memory", "Treg", "CD8 Naïve", "CD8 Memory", "CD8 Effector", "CD8 Terminally Exhausted", "γδ T lymphocytes","NK T cells", "CD56 Bright NK cells", "CD56 Dim NK cells") 


t$celltype = factor(t$celltype, levels = celltypes)
sum(is.na(t$celltype))

#t %>% 
 # group_by(ID) %>%
  #summarize(nhits = sum(p.adjust<0.05), median_NES = median(NES)) %>%
  #arrange(desc(nhits)) %>%
  #filter( grepl("KEGG", ID)) %>%
  #ggplot(aes(x=nhits)) + geom_histogram(binwidth=1) + scale_y_log10() + theme_pubr() + ggtitle("Immunologic gene sets")

#relevant_ids = t %>% 
 # group_by(ID) %>%
  #summarize(nhits = sum(p.adjust<0.05), median_NES = median(NES)) %>%
  #arrange(desc(nhits)) %>%
  #filter( grepl("KEGG", ID)) %>%
  #filter(nhits > 0) %>%
  #pull(ID) %>%
  #unique()


#p = t %>%
 # filter(ID %in% grep("KEGG", relevant_ids, value = T)) %>%
  #dplyr::select(ID, NES, p.adjust, celltype) %>%
  #mutate(
   # label = case_when(
    #  p.adjust > 0.1 ~ "",
     # p.adjust > 0.01 ~ "*",
      #p.adjust > 0.001 ~ "**",
      #!is.na(p.adjust) ~ "***",
      #TRUE ~ NA_character_
    #)) %>%
  #mutate(group = case_when(
   # celltype %in% c("CD8 Naïve",  "CD8 Effector","CD8 Memory",  "NK T cells", "γδ T lymphocytes", "CD8 Terminally Exhausted" ) ~ "CD8 T cells",
    #celltype %in% c("CD4 Naïve", "CD4 Memory", "Treg") ~ "CD4 T cells",
    #celltype %in% c("CD56 Dim NK cells",  "CD56 Bright NK cells") ~ "NK cells" )) %>%
  #ggplot(aes(y=NES, x = celltype, fill = group)) + 
  #geom_bar(stat = "identity", position = "dodge", width = 0.75, linewidth = 1.5) + 
  #facet_wrap(ID~., ncol=3) +
  #geom_text(aes(label = label, y = ifelse(NES>0, NES+0.1, NES-0.3))) + 
  #theme_pubr(base_size = 10) +
  #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.99)) + xlab("") +
  #scale_fill_manual(values = c("#d68e7e", "#f0c15d", "#9fd1cf",  "#d6bce3"), breaks = c("CD8 T cells", "CD4 T cells", "NK cells"), name = "Lineage")  
  #scale_color_manual(values = c("grey30", "brown3"))

#p

#save_plot("top_hallmark_pathways_barplot.pdf", p, base_height = 20, base_width = 30, limitsize = FALSE) 

#### Heatmap ####
relevant_ids = t %>% 
  group_by(ID) %>%
  summarize(nhits = sum(p.adjust<0.05), median_NES = median(NES)) %>%
  arrange(desc(nhits), desc(abs(median_NES))) %>%
  ungroup() %>% pull(ID) %>% unique() %>% 
  head(n=50) 

temp = t %>% 
  filter(ID %in% relevant_ids) %>%
  dplyr::select(ID, celltype, NES, p.adjust) %>%
  mutate(celltype = factor(celltype, levels = celltypes)) %>%
  arrange(celltype) %>%
  mutate(sign = case_when(p.adjust < 0.001 ~ "***",
                          p.adjust < 0.01 ~ "**",
                          p.adjust < 0.1 ~ "*",
                          .default = ""))
NES_df = temp %>%
  pivot_wider(names_from = celltype, values_from = NES, id_cols = ID)
p_val_df = temp %>%
  pivot_wider(names_from = celltype, values_from = sign, id_cols = ID)

ids = NES_df[,1]$ID
NES_df = as.matrix(NES_df[,-1])
rownames(NES_df) = ids

ids = p_val_df[,1]$ID
p_val_df = as.matrix(p_val_df[,-1])
rownames(p_val_df) = ids

mean(rownames(NES_df) == rownames(p_val_df))
mean(colnames(NES_df) == colnames(p_val_df))

# Generate factors for heatmap annotation
cohort_fac = gsub("\\+.*", "", colnames(NES_df))
celltype_fac = gsub(".*\\+", "", colnames(NES_df))
celltype_colors = c("#dddddd", rep("#f7927c", 3), rep("#ffd754", 3), rep("#8bf0ec", 1), rep("#d79df5", 3))
names(celltype_colors) = celltypes

colnames(NES_df) = celltype_fac
colnames(p_val_df) = celltype_fac

top_anno = HeatmapAnnotation(Celltype = celltype_fac,
                             annotation_name_gp = gpar(fontsize = 10),
                             border = T,
                             annotation_legend_param = list(Celltype = list(direction = "horizontal",
                                                                           at = celltypes, 
                                                                            title = "Celltypes",
                                                                            ncol =5,
                                                                            title_position = "topleft")))

col_fun = colorRamp2(c(min(NES_df), 0, max(NES_df)), c("cyan4", "white", "brown2"))

hm = Heatmap(NES_df,
             cluster_rows = F,
             cluster_columns = F,
             show_row_names = T,
             top_annotation = top_anno,
             name = "NES",
             column_names_rot = 45,
             heatmap_width = unit(28, "cm"),
             row_names_gp = gpar(fontsize = 10),
             row_title_gp = gpar(fontsize = 10),
             border = T,
             cell_fun = function(j, i, x, y, width, height, fill) {
               grid.text(sprintf("%s", p_val_df[i, j]), x, y, gp = gpar(fontsize = 8))
             })
p1 = draw(hm, heatmap_legend_side = "left", annotation_legend_side = "bottom")

# Save as a pdf
#pdf("heatmap.pdf", width = 24, height = 10)
#p1

#dev.off()
