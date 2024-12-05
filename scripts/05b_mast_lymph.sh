#!/bin/bash -l
#$ -l h_rt=96:00:00
#$ -N mast_lymph
#$ -o mast_lymph.log
#$ -m e
#$ -j y
#$ -P montilab-p
#$ -hard -l buyin=TRUE
#$ -pe omp 16
#$ -l mem_per_core=16G

cd /rprojectnb2/montilab-p/projects/brcameta/brca_atlas/scripts
module load R/4.2.1
Rscript --verbose 05b_mast_lymph.R

