#!/bin/bash -l
#$ -l h_rt=12:00:00
#$ -N imm_sccomp
#$ -o imm_sccomp_subtype.log
#$ -m e
#$ -j y
#$ -P brcameta
#$ -pe omp 16

cd /rprojectnb2/montilab-p/projects/brcameta/brca_atlas/scripts
module load R/4.2.1
Rscript --verbose 08c_imm_sccomp_subtype.R