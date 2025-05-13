# Nurefsan Sariipek
# Latest update on January 14, 2024

# GSEA based on pseudobulk DE results
# This script is based on Harvard tutorial and Ksenia's note which you can find on your local files under R/R scripts

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

# Empty environment
rm(list=ls())

# For Nurefsan:
my_wd <- "/Users/dz855/Dropbox (Partners HealthCare)/ImmuneEscapeTP53/AnalysisNurefsan/"

name = "all_celltypes"
#load the table you created end of 3_DGE script
#Make sure and pray these are pre-ranked
de_res_total = read.csv(file = paste0(my_wd, "DGE/only_Tcells-nonsig/cd8exhaustedandeffector.pseudobulk_DE_res.csv"))

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
  
  # GSE across MSigDB Hallmark immunologic gene set
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
  
# GSE across MSigDB C1 immunologic gene set
dataset = msigdbr(species = "Homo sapiens",
                  category = "C1") %>% 
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

  
 # GSE across MSigDB C2 immunologic gene set
  dataset = msigdbr(species = "Homo sapiens",
                    category = "C2") %>% 
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
  
 # GSE across MSigDB C3 immunologic gene set
  dataset = msigdbr(species = "Homo sapiens",
                    category = "C3") %>% 
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
  
  # GSE across MSigDB C4 immunologic gene set
  dataset = msigdbr(species = "Homo sapiens",
                    category = "C4") %>% 
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
  
  # GSE across MSigDB C5 immunologic gene set
  dataset = msigdbr(species = "Homo sapiens",
                    category = "C5") %>% 
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
  # GSE across MSigDB C6 immunologic gene set
  dataset = msigdbr(species = "Homo sapiens",
                    category = "C6") %>% 
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
  
# GSE across MSigDB C7 immunologic gene set
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
  
  result$celltype = i
  total_result = rbind(total_result, result)
}

write.table(total_result, paste0(name, ".total_GSEA_cd8cells_result.tsv"), col.names = T, row.names = F, quote = F, sep = '\t')

total_result %>% filter(grepl('ANTIGEN', ID)) %>% dplyr::select(-setSize, -enrichmentScore,-qvalue,-leading_edge) %>% View()

getwd()


