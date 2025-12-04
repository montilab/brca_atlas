#!/bin/bash -l
#$ -l h_rt=24:00:00
#$ -N k2_taxonomer_imm
#$ -o k2_taxonomer_imm.log
#$ -m e
#$ -j y
#$ -P brcameta
#$ -pe omp 36

cd /rprojectnb2/montilab-p/projects/brcameta/brca_atlas/scripts
module load R/4.2.1
module load libcurl
export LD_PRELOAD="${SCC_LIBCURL_LIB}/libcurl.so"
Rscript --verbose 05d_k2_imm.R
