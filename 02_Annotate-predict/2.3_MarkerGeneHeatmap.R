# Peter van Galen, 250902
# Plot a heatmap of cell type-specific marker genes

# Load libraries
library(Seurat)
library(tidyverse)
library(data.table)

# Set working directory
# fmt: skip
setwd("~/DropboxMGB/Projects/ImmuneEscapeTP53/hsct_immune_reconstitution/02_Annotate-predict")
# Or for VM
setwd("~/hsct_immune_reconstitution/02_Annotate-predict")

# Delete environment variables
rm(list = ls())

# Load data
seu <- readRDS("../AuxiliaryFiles/250528_Seurat_complete.rds")

# Define colors to use in plots
celltype_colors_df <- read.table(
  "../celltype_colors.txt",
  sep = "\t",
  header = T,
  stringsAsFactors = F,
  comment.char = ""
)
celltype_colors <- setNames(
  celltype_colors_df$color,
  celltype_colors_df$celltype
)

# Subset for cells with successful annotation
seu_subset <- subset(seu, !is.na(celltype))

# Subset for 1,000 cells per cell type so you don't need >250 GB of memory for the analyses below
set.seed(94)
cells_by <- split(colnames(seu), seu$celltype)
picked <- unlist(
  lapply(cells_by, function(v) {
    if (length(v) >= 1000) sample(v, 1000) else v
  }),
  use.names = FALSE
)
seu_subset2 <- subset(seu, cells = picked)
# Check
as_tibble(seu_subset2@meta.data) %>% count(celltype)

# Find markers for all cell types
Idents(seu_subset2) <- "celltype"
seu_subset2 <- NormalizeData(seu_subset2)
seu_subset2 <- ScaleData(seu_subset2) # Not needed for FindAllMarkers but for the heatmap
markers <- FindAllMarkers(
  seu_subset2,
  only.pos = TRUE,
  min.pct = 0.1,
  logfc.threshold = 0.1
)

# Save top 250 genes per cell type
markers_tib <- as_tibble(markers) %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 250)
write_tsv(markers_tib, file = "2.3_Celltype_markers.txt")

# Transform to a signature-like file (top 50 genes per cell type)
markers.dt.ls <- lapply(split(markers, f = markers$cluster), function(x) {
  data.table(x)
})
markers.dt.ls <- lapply(markers.dt.ls, function(x) setorder(x, -avg_log2FC))
markers.df <- do.call(cbind, lapply(markers.dt.ls, function(x) x$gene[1:50]))
write.table(
  markers.df,
  file = "2.3_Celltype_signatures.txt",
  sep = "\t",
  col.names = T,
  row.names = F,
  quote = F
)

# Manually curate markers from the signature file to visualize
# fmt: skip
heatmap_markers <- c(
"CD34", "MEIS1", "GATA2", "ITGA2B", "HDC",
"MS4A2", "TPSAB1", "TPSB2", "CPA3", "CSF2RB", "KIT",
"PF4", "PPBP", "GP9", "GP1BA", "KLF1", "HBA2", "HBB",
"GYPA", "TOP2A", "PRSS2", "SPINK2",
"CEBPA", "CSF3R", "ELANE", "PRTN3", "AZU1", "CEBPE",
"S100A8", "S100A9", "LYZ", "CD14", "CD163", "MS4A7",
"FCGR3A", "CSF1R", "CLEC10A", "CD1C", "LILRA4",
"CLEC4C", "IRF8", "RAG1", "RAG2", "VPREB1",
"IGLL1", "CD79A", "PAX5", "CD19", "MS4A1",
"IGHG1", "IGHG4", "IGHA2",
"SDC1", "CCR7", "TCF7", "LEF1", "IL7R", "IL2RA",
"FOXP3", "CTLA4", "CD8A", "CD8B",
"GZMK", "IFNG", "EOMES", "GZMH",
"NKG7", "KLRD1", "PRF1",
"IL2RB", "NCAM1", "CXCL12", "COL3A1", "COL1A1", "TF"
)
# Check
heatmap_markers %in% as.vector(markers.df)

# Create heatmap
DoHeatmap(
  seu_subset2,
  features = heatmap_markers,
  group.by = "celltype",
  group.colors = celltype_colors,
  draw.lines = F # add in Illustrator
) +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    name = "Scaled\nexpression"
  ) +
  theme(
    axis.text.y = element_text(size = 6, color = "black"),
    aspect.ratio = 0.5
  )

ggsave("2.3_MarkerGeneHeatmap.pdf", width = 20, height = 8)
