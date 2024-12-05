#!/bin/bash -l
#$ -l h_rt=8:00:00
#$ -N cellxgene
#$ -o cellxgene.log
#$ -m e
#$ -j y
#$ -P brcameta
#$ -pe omp 16

cd /rprojectnb2/montilab-p/projects/brcameta/brca_atlas/scripts
module load R/4.2.1
Rscript --verbose 09_cellxgene_prep.R