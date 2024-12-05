#!/bin/bash -l
#$ -l h_rt=72:00:00
#$ -N mast_strom
#$ -o mast_strom.log
#$ -m e
#$ -j y
#$ -P montilab-p
#$ -pe omp 16
#$ -l mem_per_core=16G

cd /rprojectnb2/montilab-p/projects/brcameta/brca_atlas/scripts
module load R/4.2.1
Rscript --verbose 05b_mast_strom.R

