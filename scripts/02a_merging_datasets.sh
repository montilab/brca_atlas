#!/bin/bash -l
#$ -l h_rt=8:00:00
#$ -N merging
#$ -o merging.log
#$ -m e
#$ -j y
#$ -P findthecause
#$ -pe omp 28

cd /rprojectnb2/montilab-p/projects/brcameta/brca_atlas/scripts
module load R/4.2.1
Rscript --verbose 02a_merging_datasets.R