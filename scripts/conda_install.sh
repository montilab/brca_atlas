#!/bin/bash -l
#$ -l h_rt=1:00:00
#$ -N conda
#$ -o conda.log
#$ -m e
#$ -j y
#$ -l gpus=1
#$ -l gpu_c=6.0
#$ -P gsc-p
#$ -pe omp 16

cd /rprojectnb2/montilab-p/projects/brcameta/brca_atlas/scripts
# For integration
mamba create -n scvi python=3.9
mamba activate scvi
mamba install -c conda-forge scanpy python-igraph leidenalg
mamba install pytorch torchvision torchaudio pytorch-cuda=12.1 -c pytorch -c nvidia
mamba install jaxlib=*=*cuda* jax cuda-nvcc -c conda-forge -c nvidia
mamba install scvi-tools -c conda-forge



# For benchmarking
mamba create -n scib python=3.9
mamba activate scib
mamba install -y scvi-tools -c conda-forge
pip install scib-metrics
pip install ipykernel
pip install --upgrade "jax[cuda12_pip]" -f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html
pip install ipykernel
python -m ipykernel install --user --name scib
