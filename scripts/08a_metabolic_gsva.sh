#!/bin/bash -l
#$ -l h_rt=6:00:00
#$ -N metabolic_gsva
#$ -o metabolic_gsva.log
#$ -m e
#$ -j y
#$ -P brcameta
#$ -pe omp 16

cd /rprojectnb2/montilab-p/projects/brcameta/brca_atlas/scripts
module load R/4.2.1
Rscript --verbose 08a_metabolic_gsva.R