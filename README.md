# TP53_ImmuneEscape

This repository contains single-cell analysis on longitudinal bone marrow samples from AML patients during the course of stem cell transplantation. We performed paired 10x Chromium Single Cell 5' and V(D)J Enrichment for 27 unique samples 44 libraries in total **<--UPDATE**.

## 0_PreProcessing
To generate count matrixes, we first created configuration files according to [10x instructions](https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/running-pipelines/cr-5p-multi) (0_PreProcessing/1_CreateConfigFiles.sh). We then submitted cellranger-8.0.1 multi jobs with standard GRCh38 and VDJ references (0_PreProcessing/2_RunCellRangerMulti.sh, 0_PreProcessing/3_Cellranger.sh).

_@nurefsan as you go through each of the steps below, please explain them in detail and with references to the original scripts_

## 1_Seurat
Load gene expression data from samples, integrate, reduce dimensionality, cluster.


## 2_Annotate
Annotate cell clusters using canonical marker genes (2.Annotate)


## 3_DGE
Add single-cell genotyping metadata to the Seurat objects containing gene expression and cell type annotations.


## 4_GSEA
Gene set enrichment analysis was performed by...

## 5_Souporcell
Distuingish cells as donor and host cells using Souporcell package (Heaton,2019)


## 6_TCR_Diversity
Analysis of TCR diversity.


## 7_Cell_Proportions



## 8_Antigen_Prediction



## 9_Numbat





