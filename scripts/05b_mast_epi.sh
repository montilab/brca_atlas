#!/bin/bash -l
#$ -l h_rt=9:00:00
#$ -N mast_epi
#$ -o mast_epi.log
#$ -m e
#$ -j y
#$ -P brcameta
#$ -t 1-10
#$ -pe omp 16
#$ -l mem_per_core=16G

cd /rprojectnb2/montilab-p/projects/brcameta/brca_atlas/scripts
module load R/4.2.1
Rscript --verbose 05b_mast_epi.R

