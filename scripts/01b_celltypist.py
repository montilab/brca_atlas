import scanpy as sc
import celltypist
from celltypist import models
import pandas as pd
import os
import re

DATA_PATH = os.path.join(os.getenv("CBM"), "otherStudies/scRNAseq/brca")
PATH = os.path.join(os.getenv("MLAB"), "projects/brcameta/brca_atlas")

task_id = int(os.environ['SGE_TASK_ID'])-1

with open("01b_h5ad_paths", 'r') as file:
    h5ad_files = [line.strip() for line in file]
    
h5ad_path = h5ad_files[task_id]

dataset_name = re.search(r'brca/([^/]+)/', h5ad_path).group(1)

print(dataset_name)

adata = sc.read_h5ad(h5ad_path)

model = models.Model.load(model = 'Cells_Adult_Breast.pkl')

predictions = celltypist.annotate(adata, model = model, majority_voting = True, mode = 'best match')

print("Done with Celltypist annotation")

predictions.predicted_labels.to_csv(os.path.join(PATH, "results/celltypist/") + dataset_name + ".csv")

adata = predictions.to_adata()

adata.write_h5ad(h5ad_path)
