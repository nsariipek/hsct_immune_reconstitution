# GSEA based on pseudobulk DE results

# Load libraries
library(SingleCellExperiment)
library(Seurat)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(viridis)
library(RColorBrewer)
library(clusterProfiler)
library(ReactomePA)
library(org.Hs.eg.db)
library(msigdbr)
library(enrichplot)

name = "all_celltypes"

de_res_total = read.csv(file = paste0(name, ".pseudobulk_DE_res.csv"))
total_result = data.frame()

for (i in unique(de_res_total$celltype)) {
  print(i)
  de_res = de_res_total %>% filter(celltype == i)
  # Sort by FC, take all genes
  gene_set = de_res %>% 
    dplyr::select(gene, log2FoldChange) %>%
    arrange(desc(log2FoldChange))
  
  eg = bitr(gene_set$gene, fromType = "SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
  # setdiff(gene_set$gene, eg$SYMBOL)
  gene_set = gene_set %>% right_join(eg, by = c("gene"="SYMBOL"))
  named_list = gene_set$log2FoldChange
  names(named_list) = gene_set$ENTREZID
  
  result = data.frame()
  
  # GSE across GO-BP categories
  dataset = msigdbr(species = "Homo sapiens",
                    category = "C5",
                    subcategory = "BP") %>% 
    dplyr::select(gs_name, entrez_gene)
  
  gsea_out = GSEA(named_list,
                  TERM2GENE = dataset,
                  pvalueCutoff = 1,
                  pAdjustMethod = "fdr",
                  verbose = FALSE,
                  eps = 0)
  if (dim(gsea_out@result)[1] > 0) {
    result = rbind(result, gsea_out@result %>% as_tibble() %>% dplyr::select(-c(core_enrichment,Description)))
  }
  
  # same for MSigDB C2:CP:KEGG, MSigDB:C7, and MSigDB:H
  
  # GSE across MSigDB C2:CP:KEGG categories
  dataset = msigdbr(species = "Homo sapiens",
                    category = "C2",
                    subcategory = "CP:KEGG") %>% 
    dplyr::select(gs_name, entrez_gene)
  
  gsea_out = GSEA(named_list,
                  TERM2GENE = dataset,
                  pvalueCutoff = 1,
                  pAdjustMethod = "fdr",
                  verbose = FALSE,
                  eps = 0)
  if (dim(gsea_out@result)[1] > 0) {
    result = rbind(result, gsea_out@result %>% as_tibble() %>% dplyr::select(-c(core_enrichment,Description)))
  }
  
  # GSE across MSigDB C2:C7 categories
  dataset = msigdbr(species = "Homo sapiens",
                    category = "C7") %>% 
    dplyr::select(gs_name, entrez_gene)
  
  gsea_out = GSEA(named_list,
                  TERM2GENE = dataset,
                  pvalueCutoff = 1,
                  pAdjustMethod = "fdr",
                  verbose = FALSE,
                  eps = 0)
  if (dim(gsea_out@result)[1] > 0) {
    result = rbind(result, gsea_out@result %>% as_tibble() %>% dplyr::select(-c(core_enrichment,Description)))
  }
  
  dataset = msigdbr(species = "Homo sapiens",
                    category = "H") %>% 
    dplyr::select(gs_name, entrez_gene)
  
  gsea_out = GSEA(named_list,
                  TERM2GENE = dataset,
                  pvalueCutoff = 1,
                  pAdjustMethod = "fdr",
                  verbose = FALSE,
                  eps = 0)
  if (dim(gsea_out@result)[1] > 0) {
    result = rbind(result, gsea_out@result %>% as_tibble() %>% dplyr::select(-c(core_enrichment,Description)))
  }
  
  dotplot(gsea_out, showCategory=20)
   gseaplot2(gsea_out, geneSetID = 1)
   
  # temp = gsea_out@result
  # temp %>% 
  #   dplyr::select(ID, NES, p.adjust) %>%
  #   filter(p.adjust<0.05) %>%
  #   arrange(desc(NES)) %>%
  #   mutate(ID = sub('REACTOME_', '', ID)) %>%
  #   mutate(ID = stringr::str_trunc(ID, 40)) %>%
  #   mutate(ID = factor(ID, levels = ID)) %>%
  #   ggplot(aes(y=NES,x=ID,fill=p.adjust)) +
  #   geom_bar(stat = "identity") +
  #   theme_pubr(base_size = 10) +
  #   theme(axis.text.x = element_text(angle=90,vjust=0.5))
  
  result$celltype = i
  total_result = rbind(total_result, result)
  
}


write.table(total_result, paste0(name, ".total_GSEA_result.tsv"), col.names = T, row.names = F, quote = F, sep = '\t')

total_result %>% filter(grepl('ANTIGEN', ID)) %>% dplyr::select(-setSize, -enrichmentScore,-qvalue,-leading_edge) %>% View()



