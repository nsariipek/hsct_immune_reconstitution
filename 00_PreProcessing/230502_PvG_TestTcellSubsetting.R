# Peter van Galen, 230503
# Check T cell subsetting and dimensionality reduction from Nurefsan's data

library(tidyverse)
library(Seurat)
library(cowplot)

# Make the objects more manageable ----------------------------------------------------------------

# Run this section only once, use diet objects next time

# First, copy data to persisent disk (in RStudio Terminal tab)
#gsutil cp gs://fc-3783b423-62ac-4c69-8c2f-98cb0ee4503b/p53_Nurefsan.rds .
#gsutil cp gs://fc-3783b423-62ac-4c69-8c2f-98cb0ee4503b/Tcellsubset.RDS .

seu <- readRDS("~/p53_Nurefsan.rds")
seu_diet <- DietSeurat(seu, dimreducs = names(seu@reductions))
saveRDS(seu_diet, file = "~/p53_Nurefsan_diet.rds")

tc <- readRDS("~/Tcellsubset.RDS")
tc_diet <- DietSeurat(tc, dimreducs = names(tc@reductions))
saveRDS(tc_diet, file = "~/Tcellsubset_diet.rds")


# Check variable genes ----------------------------------------------------------------------------

seu <- readRDS("~/p53_Nurefsan_diet.rds")
tc <- readRDS("~/Tcellsubset_diet.rds")

LabelPoints(plot = VariableFeaturePlot(seu), points = head(VariableFeatures(seu), 20), repel = T, xnudge = 0, ynudge = 0) +
  theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5)) + ggtitle("Variable genes")
LabelPoints(plot = VariableFeaturePlot(tc), points = head(VariableFeatures(tc), 20), repel = T, xnudge = 0, ynudge = 0) +
  theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5)) + ggtitle("Variable genes")

# Some plots
DimPlot(tc) + theme(aspect.ratio = 1)
p1 <- FeaturePlot(tc, features = "NCAM1") + theme(aspect.ratio = 1) + theme(aspect.ratio = 1)
p2 <- FeaturePlot(tc, features = "CD8A") + theme(aspect.ratio = 1) + theme(aspect.ratio = 1)
p3 <- FeaturePlot(tc, features = "CD4") + theme(aspect.ratio = 1) + theme(aspect.ratio = 1)
plot_grid(p1,p2,p3)

# Check variable features
VariableFeatures(tc)
tc2 <- tc
tc2 <- FindVariableFeatures(tc2, nfeatures = 2000, verbose = FALSE)
all( VariableFeatures(tc2) == VariableFeatures(tc) )


