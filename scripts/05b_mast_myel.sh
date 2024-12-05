#!/bin/bash -l
#$ -l h_rt=36:00:00
#$ -N mast_myel
#$ -o mast_myel.log
#$ -m e
#$ -j y
#$ -P apoe-signatures
#$ -pe omp 28

cd /rprojectnb2/montilab-p/projects/brcameta/brca_atlas/scripts
module load R/4.2.1
Rscript --verbose 05b_mast_myel.R

