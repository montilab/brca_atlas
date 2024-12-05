#!/bin/bash -l
#$ -l h_rt=24:00:00
#$ -N k2_taxonomer_epi
#$ -o k2_taxonomer_epi.log
#$ -m e
#$ -j y
#$ -P brcameta
#$ -pe omp 28
#$ -l mem_per_core=18G

cd /rprojectnb2/montilab-p/projects/brcameta/brca_atlas/scripts
module load R/4.2.1
Rscript --verbose 05c_k2_epi.R