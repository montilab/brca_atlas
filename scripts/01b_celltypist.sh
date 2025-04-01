#!/bin/bash -l
#$ -l h_rt=12:00:00
#$ -N celltypist_array
#$ -o celltypist.log
#$ -m e
#$ -j y
#$ -t 1-10
#$ -P brcameta
#$ -pe omp 16

cd /rprojectnb2/montilab-p/projects/brcameta/brca_atlas/scripts
conda activate celltypist
python 01b_celltypist.py