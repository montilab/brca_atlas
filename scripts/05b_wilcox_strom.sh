#!/bin/bash -l
#$ -l h_rt=6:00:00
#$ -N wilcox_strom
#$ -o wilcox_strom.log
#$ -m e
#$ -j y
#$ -P apoe-signatures
#$ -pe omp 16

cd /rprojectnb2/montilab-p/projects/brcameta/brca_atlas/scripts
module load R/4.2.1
Rscript --verbose 05b_wilcox_strom.R

