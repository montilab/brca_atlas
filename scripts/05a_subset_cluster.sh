#!/bin/bash -l
#$ -l h_rt=6:00:00
#$ -N sub_clust
#$ -o sub_clust.log
#$ -m e
#$ -j y
#$ -P apoe-signatures
#$ -pe omp 28
#$ -l mem_per_core=18G

cd /rprojectnb2/montilab-p/projects/brcameta/brca_atlas/scripts
module load R/4.2.1
Rscript --verbose 05a_subset_cluster.R