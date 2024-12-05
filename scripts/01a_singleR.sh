#!/bin/bash -l
#$ -l h_rt=24:00:00
#$ -N singleR_array
#$ -o singleR_array.log
#$ -m e
#$ -j y
#$ -t 14
#$ -P findthecause
#$ -pe omp 16
#$ -l mem_per_core=16G

cd /rprojectnb2/montilab-p/projects/brcameta/brca_atlas/scripts
module load R/4.2.1
Rscript --verbose 01a_singleR.R
