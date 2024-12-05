#!/bin/bash -l
#$ -l h_rt=24:00:00
#$ -N seurat_integration
#$ -o seurat_integration.log
#$ -m e
#$ -j y
#$ -P agedisease
#$ -pe omp 28
#$ -l mem_per_core=18G

cd /rprojectnb2/montilab-p/projects/brcameta/brca_atlas/scripts
module load R/4.2.1
Rscript --verbose 03a_seurat_integration.R
