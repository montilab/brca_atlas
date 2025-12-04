#!/bin/bash -l
#$ -l h_rt=8:00:00
#$ -N survival
#$ -o survival.log
#$ -m e
#$ -j y
#$ -P findthecause
#$ -pe omp 16

cd /rprojectnb2/montilab-p/projects/brcameta/brca_atlas/scripts
module load R/4.2.1
module load libcurl
export LD_PRELOAD="${SCC_LIBCURL_LIB}/libcurl.so"
Rscript --verbose 08d_survival_gsva.R
