
#Load the needed libraries
library(tidyverse)
library(Seurat)
library(ggplot2)
library(harmony)
library(randomcoloR)
library(RColorBrewer)
library(limma)
library(readxl)
library(cowplot)
library(SingleCellExperiment)
library(ggpubr)
library(data.table)
library(pheatmap)
library(apeglm)
library(png)
library(DESeq2)
library(Matrix)
library(ggforce)
library(cloudml)
library(Matrix.utils)

#Load the saved seurat objects

seu_diet_merged <- readRDS("/Users/dz855/Dropbox (Partners HealthCare)/ImmuneEscapeTP53/AnalysisNurefsan/GEX data/seu_diet_merged.rds")

View(seu_diet_merged@meta.data)

ex <- subset(x= seu_diet_merged, subset =cohort %in% c("cohort1","cohort2"))
ex <- subset(x= ex, subset =status %in% c("remission"))
ex <- subset(x=ex, subset = id %in% c("P01.1Rem", "P01.1RemT", "P01.2Rem", "P02.1Rem", "P04.1Rem", "P04.1RemT", "P05.1Rem", "P06.1Rem", "P07.1Rem", "P07.1RemT", "P08.1Rem", "P08.1RemT"))

#Pool MNC and CD3 samples
#ex <- ex%>% mutate(pt_timepoint = gsub("RemT", "Rem", id)) 

#remove the . for the pseudobulk analysis
ex$id =gsub("\\.", "", ex$id)
#ex$groups = NULL
View(ex@meta.data)
sce = as.SingleCellExperiment(ex)
dim(counts(sce))
counts(sce)[1:6, 1:6]

dim(colData(sce))
head(colData(sce))

celltype_names = unique(sce$celltype)

sce$sample_id = paste0(sce$id)

#Subset metadata to include only the variables you want to aggregate across (here, we want to aggregate by sample and by celltype)

groups = colData(sce)[, c("celltype", "sample_id")]
head(groups)

# Aggregate across cluster-sample groups
# transposing row/columns to have cell_ids as row names matching those of groups
aggr_counts = aggregate.Matrix(t(counts(sce)), 
                               groupings = groups, fun = "sum") 


# Explore output matrix
class(aggr_counts)
dim(aggr_counts)
aggr_counts[1:6, 1:6]

# Transpose aggregated matrix to have genes as rows and samples as columns
aggr_counts = t(aggr_counts)
aggr_counts[1:6, 1:6]

tstrsplit(colnames(aggr_counts), "_") %>% str()


head(colnames(aggr_counts), n = 10)
head(tstrsplit(colnames(aggr_counts), "_")[[1]], n = 10)

celltype_names
## Initiate empty list
counts_ls = list()

#put in a loop
for (i in 1:length(celltype_names)) {
  
  ## Extract indexes of columns in the global matrix that match a given cluster
  column_idx = which(tstrsplit(colnames(aggr_counts), "_")[[1]] == celltype_names[i])
  
  ## Store corresponding sub-matrix as one element of a list
  counts_ls[[i]] = aggr_counts[, column_idx]
  names(counts_ls)[i] = celltype_names[i]
  
}
str(counts_ls)

# Explore the different components of the list
str(counts_ls)

head(colData(sce))

View(ex@meta.data)

# Extract sample-level variables
metadata = colData(sce) %>% 
  as.data.frame() %>% 
  dplyr::select(patient_identity,cohort, status, sample_id)

metadata = metadata[!duplicated(metadata), ]

dim(metadata)
head(metadata)

rownames(metadata) = metadata$sample_id
head(metadata)

t <- table(colData(sce)$sample_id,
           colData(sce)$celltype)
t[1:6, 1:6]

metadata_ls <- list()

for (i in 1:length(counts_ls)) {
  
  ## Initiate a data frame for cluster i with one row per sample (matching column names in the counts matrix)
  df <- data.frame(celltype_sample_id = colnames(counts_ls[[i]]))
  
  ## Use tstrsplit() to separate cluster (cell type) and sample IDs
  df$celltype <- tstrsplit(df$celltype_sample_id, "_")[[1]]
  df$sample_id  <- tstrsplit(df$celltype_sample_id, "_")[[2]]
  
  
  ## Retrieve cell count information for this cluster from global cell count table
  idx <- which(colnames(t) == unique(df$celltype))
  cell_counts <- t[, idx]
  
  ## Remove samples with zero cell contributing to the cluster
  cell_counts <- cell_counts[cell_counts > 0]
  
  ## Match order of cell_counts and sample_ids
  sample_order <- match(df$sample_id, names(cell_counts))
  cell_counts <- cell_counts[sample_order]
  
  ## Append cell_counts to data frame
  df$cell_count <- cell_counts
  
  ## Join data frame (capturing metadata specific to cluster) to generic metadata
  df <- plyr::join(df, metadata, 
                   by = intersect(names(df), names(metadata)))
  
  ## Update rownames of metadata to match colnames of count matrix, as needed later for DE
  rownames(df) <- df$cluster_sample_id
  
  ## Store complete metadata for cluster i in list
  metadata_ls[[i]] <- df
  names(metadata_ls)[i] <- unique(df$celltype)
  
}

# Explore the different components of the list
str(metadata_ls)


#### DGE ####
# Select cell type of interest
celltype_names

all(names(counts_ls) == names(metadata_ls))
##
#Run this when you are running one celltype ----- for a test run#
idx <- which(names(counts_ls) == "CD4 Memory")
cluster_counts <- counts_ls[[idx]]
cluster_metadata <- metadata_ls[[idx]]

# Check contents of extracted objects
cluster_counts[1:6, 1:6]
head(cluster_metadata)

# Check matching of matrix columns and metadata rows
all(colnames(cluster_counts) == rownames(cluster_metadata))

#run this code to correct rownames when running for only one cell type
cluster_metadata = cluster_metadata %>% `row.names<-`(cluster_metadata$celltype_sample_id)

#----------#
#Ksenia's code- for the loop of all the celltypes
total_res_table_thres = data.frame()

for (celltype in unique(sce$celltype)) {
  print(celltype)
  
  # Only run DE if cell counts for both conditions exceed 99
  cell_counts = metadata_ls[[celltype]] %>% as_tibble() %>% group_by(cohort) %>% summarize(n=sum(cell_count)) %>% mutate(ind = n >=100) %>% pull(ind) %>% sum()
  if (cell_counts < 2) {
    print(paste0(celltype, ' doesn\'t have enough cells, skip'))
    next
  }
  idx = which(names(counts_ls) == celltype)
  cluster_counts = counts_ls[[idx]]
  cluster_metadata = metadata_ls[[idx]]
# Create DESeq2 object 
dds = DESeqDataSetFromMatrix(cluster_counts, 
                             colData = cluster_metadata, 
                             design =  ~ cohort)

#new line from Ksenia to make sure that the plots are generated according to correct reference
ddsColl$cohort = relevel(dds$cohort, ref = "cohort1")
# Transform counts for data visualization
rld = rlog(dds, blind=TRUE)

# Plot PCA
p1 = DESeq2::plotPCA(rld, ntop = 500, intgroup = "cohort") + ggrepel::geom_text_repel(aes(label = name)) + theme_pubr()
p1

DESeq2::plotPCA(rld, ntop = 500, intgroup = "cell_count")

# Extract the rlog matrix from the object and compute pairwise correlation values
rld_mat = assay(rld)
rld_cor = cor(rld_mat)


annotation = cluster_metadata[, c("cohort"), drop=F] %>% rownames_to_column("temp") %>% mutate(temp = sub("-.\\.", ".", temp)) %>% unique() 
rownames(annotation) = annotation$temp


# Run DESeq2 differential expression analysis
dds = DESeq(dds, quiet = T)


# Generate results object
res = results(dds, 
               name = "cohort_cohort2_vs_cohort1",
               alpha = 0.05)

# Shrink the log2 fold changes to be more appropriate using the apeglm method - should cite [paper]() when using this method
res = lfcShrink(dds, 
                coef = "cohort_cohort2_vs_cohort1",
                res=res,
                type = "apeglm",
                quiet = T)

res_tbl <- res %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  arrange(padj)


# Set thresholds
padj_cutoff = 0.05

# Subset the significant results
sig_res = dplyr::filter(res_tbl, padj < padj_cutoff) %>%
  dplyr::arrange(padj)

if (nrow(sig_res) < 2) {
  print(paste0(celltype, ' doesn\'t have enough DEGs, skip'))
  next
}

## Extract normalized counts from ddsColl object
normalized_counts = counts(dds, normalized = TRUE)

## Extract top 20 DEG from resLFC (make sure to order by padj)
top20_sig_genes = sig_res %>%
  dplyr::arrange(padj) %>%
  dplyr::pull(gene) %>%
  head(n = 20)

## Extract matching normalized count values from matrix
top20_sig_counts = normalized_counts[rownames(normalized_counts) %in% top20_sig_genes, ]

## Convert wide matrix to long data frame for ggplot2
top20_sig_df = data.frame(top20_sig_counts)
top20_sig_df$gene = rownames(top20_sig_counts)

top20_sig_df = melt(setDT(top20_sig_df), 
                    id.vars = c("gene"),
                    variable.name = "sample") %>% 
  data.frame()

## Replace "." by " " in cluster_sample_id variable (melt() introduced the ".")
top20_sig_df$sample= gsub("\\.", " ", top20_sig_df$sample)


top20_sig_df$celltype_sample_id= top20_sig_df$sample

top20_sig_df = plyr::join(top20_sig_df, as.data.frame(colData(dds)),
                          by = "celltype_sample_id")

## Generate plot
p3 = ggplot(top20_sig_df, aes(y = value, x = cohort, col = cohort)) +
  geom_jitter(height = 0, width = 0.15) +
  scale_y_continuous(trans = 'log10') +
  ylab("log10 of normalized expression level") +
  xlab("condition") +
  ggtitle("Top 20 Significant DE Genes") +
  theme(plot.title = element_text(hjust = 0.5)) +
  facet_wrap(~ gene)


###I don't know wwhy but heatmap crashes the code so I removed the heatmap part. 

## Extract normalized counts for significant genes only
sig_counts = normalized_counts[rownames(normalized_counts) %in% sig_res$gene, ]

## Set a color-blind friendly palette
heat_colors = rev(brewer.pal(11, "PuOr"))


## Run pheatmap using the metadata data frame for the annotation

# Volcano plot
log2fc_cutoff = 0.26 # (20% increase)
res_table_thres = res_tbl[!is.na(res_tbl$padj), ] %>% 
  mutate(threshold = padj < padj_cutoff & abs(log2FoldChange) >= log2fc_cutoff)

top_genes = top20_sig_df$gene %>% unique()

## Generate plot
p5 = ggplot() +
  geom_point(data=res_table_thres, aes(x = log2FoldChange, y = -log10(padj), colour = threshold)) +
  ggrepel::geom_text_repel(data=res_table_thres %>% filter(gene %in% top_genes), aes(x = log2FoldChange, y = -log10(padj), label = gene)) +
  xlab("log2 fold change") +
  ylab("-log10 adjusted p-value") +
  scale_color_manual(values = c("grey60", "red3")) +
  theme_pubr() + 
  ggtitle("DE")

res_table_thres$celltype = celltype
total_res_table_thres = rbind(total_res_table_thres, res_table_thres)
p = ggarrange(ggarrange(p1,p5, ncol=1, widths = c(0.3,0.7), align = "none"),p3, align = "none")
save_plot(paste0("~/", sub('/', '_', celltype), ".summary.pdf"), p, base_height = 20, base_width = 20)

}
​
### Write all results to file
name = "all_celltypes"
write.csv(x = total_res_table_thres, file = paste0(name, ".pseudobulk_DE_res.csv"), quote = F, row.names = F)






