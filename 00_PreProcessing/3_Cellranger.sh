#!/bin/bash
#$ -l h_vmem=16G
#$ -pe smp 4
#$ -binding linear:4
#$ -l h_rt=12:00:00
#$ -cwd
# This is required to use dotkits inside scripts
source /broad/software/scripts/useuse
use .cellranger-8.0.1

# Assign the first argument to the sample_id variable
sample_id="$1"
seqdate="$2"

# Start run
cellranger multi --id="$sample_id" --csv="/broad/vangalenlab/sariipek/cellranger/multi_configs_${seqdate}/${sample_id}_config.csv"

