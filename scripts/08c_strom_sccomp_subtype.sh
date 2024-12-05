#!/bin/bash -l
#$ -l h_rt=12:00:00
#$ -N strom_sccomp
#$ -o strom_sccomp_subtype.log
#$ -m e
#$ -j y
#$ -P brcameta
#$ -pe omp 16

cd /rprojectnb2/montilab-p/projects/brcameta/brca_atlas/scripts
module load R/4.2.1
Rscript --verbose 08c_strom_sccomp_subtype.R