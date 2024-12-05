#!/bin/bash -l
#$ -l h_rt=18:00:00
#$ -N scib-metrics_combined
#$ -o scib-metrics_combined.log
#$ -m e
#$ -j y
#$ -P montilab-p
#$ -l gpus=1
#$ -l gpu_c=6.0
#$ -pe omp 16

cd /rprojectnb2/montilab-p/projects/brcameta/brca_atlas/scripts
conda activate scib
python 04_scib_benchmark.py