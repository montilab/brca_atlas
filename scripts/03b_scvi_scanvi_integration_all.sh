#!/bin/bash -l
#$ -l h_rt=48:00:00
#$ -N scvi_scanvi
#$ -o scvi_scanvi_all.log
#$ -m e
#$ -l gpus=1
#$ -l gpu_c=6.0
#$ -j y
#$ -P yaptaz
#$ -pe omp 16

cd /rprojectnb2/montilab-p/projects/brcameta/brca_atlas/scripts
conda activate scvi
python 03b_scvi_scanvi_integration_all.py