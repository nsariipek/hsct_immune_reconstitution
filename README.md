# TP53_ImmuneEscape

This repository contains single-cell analysis from 33 AML patients during the course of stem cell transplantation. We performed paired 10x Chromium Single Cell 5' and V(D)J Enrichment for 48 unique samples and 65 libraries in total. Of these, 12 patients had longitudinal samples; meanwhile, 21 patients only had 100 days post transplant samples. 

## 0_PreProcessing
To generate count matrixes, we first created configuration files according to [10x instructions](https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/running-pipelines/cr-5p-multi) (0_PreProcessing/1_CreateConfigFiles.sh). We then submitted cellranger-8.0.1 multi jobs with standard GRCh38 and VDJ references (0_PreProcessing/2_RunCellRangerMulti.sh, 0_PreProcessing/3_Cellranger.sh).

## 1_Seurat
Load and integrate gene expression data from multiple samples, ensuring proper sample identification. Assess cell numbers, visualize quality control (QC) metrics, and apply appropriate filtering. (1.1_CreateSeuratObject.R)
(1.2_DimensionalityReduction_Clustering)Given the large dataset, we split the Seurat object into 20% and 80% subsets. Normalization and dimensionality reduction were performed on the 20% subset. We then constructed the nearest neighbor graph, performed clustering, and generated the initial UMAP visualization. Annotation was conducted using only the 20% subset of the data. (Refer to 1.2_DimensionalityReduction_Clustering for details.)

## 2_Annotate
Load the 20% subset of the Seurat object and annotate cell clusters using canonical marker genes. Additionally, we utilized three distinct gene sets from our lab—two unpublished and one published—to compute module scores. These scores were visualized using feature plots to help with annotations. (2.1_Annotate.R)
Subsetted T and NK cells for a more granular analysis, re-identified variable features, performed clustering, and generated a UMAP specifically for the T cell subset.For T cell annotations, we utilized one unpublished annotation from our lab, along with several published studies and curated gene sets for each T and NK cell subset. (2.2_Tcell_Annotation.R)
After completing the annotations, we merged the 20% and 80% subsets of the Seurat object, using the annotated 20% section as a reference. This process followed the tutorial from Seurat's integration and mapping guide. Finally, we visualized the entire dataset using UMAP.(2.4_Merge_497K_cells.R)

## 3_DGE
Add single-cell genotyping metadata to the Seurat objects containing gene expression and cell type annotations.


## 4_GSEA
Gene set enrichment analysis was performed by...

## 5_Souporcell
Distinguish cells as donor and host cells using the Souporcell package (Heaton,2019)
First, we extracted the barcodes from the Seurat object for each patient; for longitudinal patients, we merged the barcodes and removed the duplicates across different time points. (5.1_Extracting Barcodes for the Souporcell) 
In the second part, using samtools, we merged the bam files for the patients who had multiple longitudinal samples. For the samples that did not have the longitudinal samples, we just used the output processed bam file.(5.2_MergeBamfiles)


## 6_TCR_Diversity
Analysis of TCR diversity.
Begin by loading the cell ranger multioutput /filtered_contig_annotations.csv for each sample and start combining and integrating them with T cell gene expression data. We used the scRepertoire(https://github.com/BorchLab/scRepertoire) package for this analysis but modified the diversity calculation approach within the script. (6.1_PostTx_3to6months.R) In this script, we focused specifically on 100-day post-treatment samples, comparing groups with different outcomes using parameters derived from single-cell expression data. The output is saved for further analysis. Additionally, you can subset the data based on cell type annotations or Souporcell output for more distinct visualizations.


## 7_Cell_Proportions
For longitudinal samples, we visualized cell proportions using alluvial plots to illustrate their dynamic changes over time (7.1_AlluvialPlots.R).
To analyze the CD4/CD8 ratio and markers indicative of a healthy environment, we used a separate script (7.2_CD4/8ratio.R).
Finally, by incorporating Souporcell results, we visualized recipient and donor cell ratios across different time points (7.3_Donor/Host.R).


## 8_Antigen_Prediction



## 9_Numbat





