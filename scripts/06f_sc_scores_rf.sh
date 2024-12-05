#!/bin/bash -l
#$ -l h_rt=24:00:00
#$ -N epi_scoring_rf
#$ -o epi_scoring_rf.log
#$ -m e
#$ -j y
#$ -P brcameta
#$ -pe omp 16

cd /rprojectnb2/montilab-p/projects/brcameta/brca_atlas/scripts
module load R/4.2.1
Rscript --verbose 06f_sc_scores_rf.R