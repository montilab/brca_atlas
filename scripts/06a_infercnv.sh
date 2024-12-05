#!/bin/bash -l
#$ -l h_rt=48:00:00
#$ -N infercnv
#$ -o infercnv.log
#$ -m e
#$ -j y
#$ -t 1-8
#$ -P montilab-p
#$ -pe omp 28

# Define the specific task IDs you want to run
# desired_tasks=(4 7)

# Check if the current task ID is in the desired list
# if [[ " ${desired_tasks[@]} " =~ " ${SGE_TASK_ID} " ]]; then
#   # Execute your task here
#   echo "Running task ID: ${SGE_TASK_ID}"
#   # Add your task-specific commands here
# else
#   echo "Skipping task ID: ${SGE_TASK_ID}"
# fi

cd /rprojectnb2/montilab-p/projects/brcameta/brca_atlas/scripts
module load R/4.2.1
Rscript --verbose 06a_infercnv.R