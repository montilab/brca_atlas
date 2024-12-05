#!/bin/bash -l
#$ -l h_rt=12:00:00
#$ -N k2_taxonomer_strom
#$ -o k2_taxonomer_strom.log
#$ -m e
#$ -j y
#$ -P brcameta
#$ -pe omp 16

cd /rprojectnb2/montilab-p/projects/brcameta/brca_atlas/scripts
module load R/4.2.1
Rscript --verbose 05c_k2_strom.R