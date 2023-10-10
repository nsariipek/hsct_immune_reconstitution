#Load the saved TP53 data set. This already has undergone some standard QC filtering steps.
seu_diet_merged <- readRDS("~/seu_diet_merged.rds")
tc_diet<- readRDS("~/Tcellsubset_diet.rds")

#check the metadata
View(seu_diet_merged@meta.data)  
View(tc_diet@meta.data)

#This dataset includes all of the samples. If you want to subset a certain time period, do it after this step.
#subset only 0-6 months post-tx samples from cohorts 1-2
tcremonly <-  subset(x= tc_diet, subset =cohort %in% c("cohort1","cohort2"))
tcremonly <- subset(x= tcremonly, subset =status %in% c("remission"))
tcremonly2_6mo <- subset(x=tcremonly, subset = id %in% c("P01.1Rem", "P01.1RemT", "P01.2Rem", "P02.1Rem", "P04.1Rem", "P04.1RemT", "P05.1Rem", "P06.1Rem", "P07.1Rem", "P07.1RemT", "P08.1Rem", "P08.1RemT"))



