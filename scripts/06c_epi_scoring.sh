#!/bin/bash -l
#$ -l h_rt=24:00:00
#$ -N epi_scoring
#$ -o epi_soring.log
#$ -m e
#$ -j y
#$ -P brcameta
#$ -pe omp 16

cd /rprojectnb2/montilab-p/projects/brcameta/brca_atlas/scripts
module load R/4.2.1
Rscript --verbose 06c_epi_scoring.R