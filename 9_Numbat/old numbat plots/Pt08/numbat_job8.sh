#!/bin/bash
#$ -l h_vmem=64G
#$ -pe smp 1
#$ -binding linear:1
#$ -l h_rt=30:00:00
#$ -cwd

source /broad/software/scripts/useuse
use R-4.1

Rscript Numbat_script.R
