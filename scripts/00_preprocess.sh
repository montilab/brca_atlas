#!/bin/bash -l
#$ -l h_rt=4:00:00
#$ -N preprocess
#$ -o preprocess.log
#$ -m e
#$ -j y
#$ -P gsc-p
#$ -pe omp 16

cd /rprojectnb2/montilab-p/projects/brcameta/brca_atlas/scripts
module load R/4.2.1
Rscript --verbose 00_preprocess.R