#!/bin/bash -l
#$ -l h_rt=12:00:00
#$ -N adding_embeddings
#$ -o adding_embeddings.log
#$ -m e
#$ -j y
#$ -P agedisease
#$ -pe omp 28
#$ -l mem_per_core=18G

cd /rprojectnb2/montilab-p/projects/brcameta/brca_atlas/scripts
module load R/4.2.1
Rscript --verbose 03c_adding_embeddings.R