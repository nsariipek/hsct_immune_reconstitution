# Nurefsan attemps to try TREX
# 240314
# Load libraries

library(Trex)
library(ggplot2)
library(tidyr)
library(Seurat)
library(viridis)
library(tensorflow)
library(reticulate)
library(janitor)

# Clean the environment
rm(list=ls())

# Load your Seurat object:
my_wd <- "/Users/dz855/Dropbox (Partners HealthCare)/ImmuneEscapeTP53/"

#Visualize the UMAP
DimPlot(Tcells_combined, reduction = "umap", repel = T, group.by = "celltype", label = T) + theme(aspect.ratio = 1)

#maTREX will use the combined outpit list from the combineTCR()
#runTREX needs a single cell object
Tx <- runTrex(Tcells_combined, 
                     chains = "TRB",
                     encoder.model = "VAE", 
                     encoder.input = "KF",
                     reduction.name = "Trex.KF")

#Generating UMAP from Trex Neighbors
Tx <- RunUMAP(Tx, 
                     reduction = "Trex.KF",
                     dims = 1:30,
                     reduction.name = 'Trex.umap', 
                     reduction.key = 'trexUMAP_')



#Trex UMAP
plot1 <- DimPlot(Tx, 
                 reduction = "Trex.umap") + NoLegend()
plot2 <- DimPlot(Tx, 
                 group.by = "CTaa", 
                 reduction = "Trex.umap") + 
                 scale_color_viridis(discrete = TRUE, option = "B")  + 
                 theme(plot.title = element_blank()) +
                 NoLegend()

plot1 + plot2

#In order to generate a single-cell object based on the CoNGA approach, Trex offers the function CoNGAfy(). For method, select either “mean” or “dist” as described above. After performing CoNGAfy(), the user can use any of the above reduction strategies.
library(dplyr, include.only = c("%>%"))

CoNGA.seurat <- CoNGAfy(Tx, 
                        method = "dist")

CoNGA.seurat <- runTrex(CoNGA.seurat, 
                        chains = "TRB",
                        encoder.model = "VAE", 
                        encoder.input = "KF",
                        reduction.name = "Trex.KF")

CoNGA.seurat <- CoNGA.seurat %>%
  FindNeighbors(reduction = "umap") %>%
  FindClusters(algorithm = 3)

CoNGA.seurat <- RunUMAP(CoNGA.seurat, 
                        reduction = "Trex.KF", 
                        dims = 1:20, 
                        reduction.name = 'Trex.umap', 
                        reduction.key = 'trexUMAP_')

DimPlot(CoNGA.seurat, reduction = "Trex.umap") + NoLegend()


#Towards find epitope relationships, Trex has a built in data base of TCRA and TCRB sequences associated with epitopes. To append the database to the single-cell object (either before or after CoNGAfy()), you can use annotateDB()
#CoNGA.seurat <- annotateDB(CoNGA.seurat, 
                           #chains = "TRB")

#DimPlot(CoNGA.seurat, 
        #reduction = "umap", 
        #group.by = "TRB_Epitope.species") + 
        #theme(plot.title = element_blank())


#DimPlot(CoNGA.seurat, 
       # reduction = "umap",
       # group.by = "topspecies") +
  #theme(plot.title = element_blank(), aspect.ratio = 1)

#Or annotateDB() can be used on the full single-cell object and examine the sequence information along the RNA-based UMAP. An added feature to the function allows the annotations to be extended to CDR3 sequences that are within defined edit distances from the reference using edit.distance.

Tx <- annotateDB(Tx, chains = "TRB", 
                     edit.distance = 2)

top <- c("CMV","Immunization","M.tuberculosis", "Influenza", "SelfAg", "InfluenzaA", "EBV", "CMV;Immunization", "SARS-CoV-2", "CMV;EBV", "HIV")

Tx@meta.data$topspecies <- ifelse(Tx@meta.data$TRB_Epitope.species %in% top, paste("",Tx@meta.data$TRB_Epitope.species), "Others")

DimPlot(Tx, 
        reduction = "umap", 
        group.by = "topspecies", split.by = "cohort", 
        col= color_values) + 
  theme(plot.title = element_blank())
# save seurat object to use another time

saveRDS(Tx,"antigen.prediction.rds")
