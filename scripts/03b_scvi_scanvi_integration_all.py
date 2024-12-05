import scanpy as sc
import scanpy.external as sce
import scvi
import torch
from scipy import sparse
import numpy as np
import pandas as pd
from rich import print
import os

torch.set_float32_matmul_precision("high")
scvi.settings.seed = 0
print("Last run with scvi-tools version:", scvi.__version__)

DATA_PATH = os.path.join(os.getenv("CBM"), "otherStudies/scRNAseq/brca")
PATH = os.path.join(os.getenv("MLAB"), "projects/brcameta/brca_atlas")

# adata = sc.read_h5ad(os.path.join(PATH, "data/sc/combined_anndata_pp.h5ad"))
adata = sc.read_h5ad(os.path.join(PATH, "data/sc/combined_anndata.h5ad"))

# # Adding subtype + cell type category
# adata.obs['subtype_celltype'] = adata.obs['subtype_new'].astype(str) + "_" + adata.obs['celltypist_pred'].astype(str)
# print(adata.obs['subtype_celltype'].value_counts())

# # Data Processing
# adata.raw = adata  # keep full dimension safe
# adata.layers["counts"] = adata.X # Add raw counts as a layer
# adata.X = sparse.csr_matrix(adata.layers["counts"])

# # Using variable genes from Seurat
# assert(sum(~adata.var['var.features'].isna()) == 5000)
# adata.var = adata.var[['var.features']]
# adata.var['highly_variable'] = ~adata.var['var.features'].isna()

# # No integration
# sc.tl.pca(adata, n_comps=200, use_highly_variable=True)
# adata.write_h5ad(os.path.join(PATH, "data/sc/combined_anndata.h5ad"))

adata = adata[:, adata.var.highly_variable].copy()
adata.obsm["Unintegrated"] = adata.obsm["X_pca"]

# Integration with scvi
# print("Starting scvi")
# scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key="batch")
# vae = scvi.model.SCVI(adata, n_layers=2, n_latent=30, gene_likelihood="nb")
# vae.train()
# vae.save(os.path.join(PATH, "data/models/all/scvi_30"), overwrite=True, save_anndata=True)
# vae = scvi.model.SCVI.load(os.path.join(PATH, "data/models/all/scvi_30"), adata)
# SCVI_LATENT_KEY = "X_scVI"
# adata.obsm[SCVI_LATENT_KEY] = vae.get_latent_representation()
# np.savetxt(os.path.join(PATH, "data/embeddings/all/scvi30.csv"), adata.obsm[SCVI_LATENT_KEY], delimiter=",")

# Integration with scanvi (subtype x celltypist)
print("Starting subtype x celltypist Scanvi")
vae = scvi.model.SCVI.load(os.path.join(PATH, "data/models/all/scvi_30"), adata)
lvae_subtype = scvi.model.SCANVI.from_scvi_model(
    vae,
    adata=adata,
    labels_key="subtype_celltype",
    unlabeled_category="Unassigned",
)
lvae_subtype.train(max_epochs=20, n_samples_per_label=100)
lvae_subtype.save(os.path.join(PATH, "data/models/all/scanvi_subtype_celltypist_30"), overwrite=True, save_anndata=True)
SCANVI_LATENT_KEY = "X_scANVI_subtype_broad"
adata.obsm[SCANVI_LATENT_KEY] = lvae_subtype.get_latent_representation(adata)
np.savetxt(os.path.join(PATH, "data/embeddings/all/scanvi_30_subtype_celltypist.csv"), adata.obsm[SCANVI_LATENT_KEY], delimiter=",")

# Integration with scanvi (celltypist)
print("Starting celltypist Scanvi")
vae = scvi.model.SCVI.load(os.path.join(PATH, "data/models/all/scvi_30"), adata)
adata.obs['celltypist_pred'] = adata.obs['celltypist_pred'].astype('str')
adata.obs['celltypist_pred'] = adata.obs['celltypist_pred'].fillna("Unassigned")
print(adata.obs['celltypist_pred'].value_counts())

lvae_subtype = scvi.model.SCANVI.from_scvi_model(
    vae,
    adata=adata,
    labels_key="celltypist_pred",
    unlabeled_category="Unassigned",
)
lvae_subtype.train(max_epochs=20, n_samples_per_label=100)
lvae_subtype.save(os.path.join(PATH, "data/models/all/scanvi_celltypist_30"), overwrite=True, save_anndata=True)
SCANVI_LATENT_KEY = "X_scANVI_celltypist"
adata.obsm[SCANVI_LATENT_KEY] = lvae_subtype.get_latent_representation(adata)
np.savetxt(os.path.join(PATH, "data/embeddings/all/scanvi_30_celltypist.csv"), adata.obsm[SCANVI_LATENT_KEY], delimiter=",")

# Integration with Scanorama
# print("Starting Scanorama")
# idx = adata.obs.sort_values("batch").index
# adata_new = adata[idx,]
# adata_new.X = adata_new.X.toarray()
# sce.pp.scanorama_integrate(adata_new, 'batch', batch_size=2500)
# adata.write_h5ad(os.path.join(PATH, "data/sc/combined_anndata.h5ad"))
# print('X_scanorama' in adata.obsm)
# adata.obsm["Scanorama_MDE"] = scvi.model.utils.mde(adata.obsm["X_scanorama"])
