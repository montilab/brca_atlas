import scanpy as sc
import scvi
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from rich import print
from scib_metrics.benchmark import Benchmarker, BioConservation, BatchCorrection
from scvi.model.utils import mde
import os
import jax
import pickle

os.environ['XLA_PYTHON_CLIENT_PREALLOCATE'] ='false'
os.environ['XLA_PYTHON_CLIENT_ALLOCATOR']='platform'
os.environ['TF_FORCE_GPU_ALLOW_GROWTH'] = 'true'

print(jax.devices())

DATA_PATH = os.path.join(os.getenv("CBM"), "otherStudies/scRNAseq/brca")
PATH = os.path.join(os.getenv("MLAB"), "projects/brcameta/brca_atlas")

adata = sc.read_h5ad(os.path.join(PATH, "data/sc/combined_anndata.h5ad"))

# Loading Embeddings
harmony_embeddings = pd.read_csv(os.path.join(PATH, "data/embeddings/all/harmony_embedding_combined.csv"), index_col=0)
rpca_embeddings = pd.read_csv(os.path.join(PATH, "data/embeddings/all/rpca_embedding_combined.csv"), index_col=0)
mnn_embeddings = pd.read_csv(os.path.join(PATH, "data/embeddings/all/mnn_embedding_combined.csv"), index_col=0)
scvi = pd.read_csv(os.path.join(PATH, "data/embeddings/all/scvi.csv"), header=None)
scvi_30 = pd.read_csv(os.path.join(PATH, "data/embeddings/all/scvi30.csv"), header=None)
scanvi_ct = pd.read_csv(os.path.join(PATH, "data/embeddings/all/scanvi_celltypist.csv"), header=None)
scanvi_30_ct = pd.read_csv(os.path.join(PATH, "data/embeddings/all/scanvi_30_celltypist.csv"), header=None)
scanvi_subtype_celltypist = pd.read_csv(os.path.join(PATH, "data/embeddings/all/scanvi_subtype_celltypist.csv"), header=None)
scanvi_30_subtype_celltypist = pd.read_csv(os.path.join(PATH, "data/embeddings/all/scanvi_30_subtype_celltypist.csv"), header=None)

np.array_equal(harmony_embeddings.index, mnn_embeddings.index)
np.array_equal(harmony_embeddings.index, rpca_embeddings.index)
np.array_equal(harmony_embeddings.index, adata.obs.index)

adata.obsm['X_seurat_harmony'] = harmony_embeddings
adata.obsm['X_seurat_rpca'] = rpca_embeddings
adata.obsm['X_seurat_mnn'] = mnn_embeddings
adata.obsm["X_scVI"] = scvi.values
adata.obsm["X_scVI_30"] = scvi_30.values
adata.obsm["X_scANVI_ct"] = scanvi_ct.values
adata.obsm["X_scANVI_30_ct"] = scanvi_30_ct.values
adata.obsm["X_scANVI_sub_ct"] = scanvi_subtype_celltypist.values
adata.obsm["X_scANVI_30_sub_ct"] = scanvi_30_subtype_celltypist.values
adata.obsm['X_concat'] = np.concatenate((adata.obsm['X_seurat_harmony'], 
                                         adata.obsm['X_seurat_rpca'],
                                         adata.obsm['X_scVI_30'],
                                         adata.obsm["X_scANVI_30_ct"]), axis = 1)
adata.obsm["Unintegrated"] = adata.obsm["X_pca"]

# Finding unique embeddings
unique_idx = []
for embed in ["Unintegrated", "X_scANVI_ct", "X_scANVI_30_ct", "X_scANVI_sub_ct", 
"X_scANVI_30_sub_ct", "X_scVI", "X_scVI_30", "X_seurat_harmony", 
"X_seurat_rpca", "X_seurat_mnn", 'X_concat']:
    print(embed)
    print(adata.obsm[embed].shape)
    print(np.unique(adata.obsm[embed], axis=0).shape)
    unique_idx.append(np.unique(adata.obsm[embed], axis=0, return_index=True)[1])

from functools import reduce
unique_idx = reduce(np.intersect1d, unique_idx)
adata_sub = adata[unique_idx, :]

# Evaluation
biocons = BioConservation(isolated_labels=True, nmi_ari_cluster_labels_leiden=False, nmi_ari_cluster_labels_kmeans=True, silhouette_label=True, clisi_knn=True)
# biocons = BioConservation(isolated_labels=False, nmi_ari_cluster_labels_leiden=False, nmi_ari_cluster_labels_kmeans=False, silhouette_label=False, clisi_knn=False)
batchcons = BatchCorrection(silhouette_batch=True, ilisi_knn=True, kbet_per_label=True, graph_connectivity=True, pcr_comparison=False)
bm = Benchmarker(
    adata_sub,
    batch_key="batch",
    label_key="celltypist_pred",
    bio_conservation_metrics=biocons,
    batch_correction_metrics=batchcons,
    embedding_obsm_keys=["Unintegrated", "X_scANVI_ct", "X_scANVI_30_ct", "X_scANVI_sub_ct", 
"X_scANVI_30_sub_ct", "X_scVI", "X_scVI_30", "X_seurat_harmony", 
"X_seurat_rpca", "X_seurat_mnn", "X_concat"],
    pre_integrated_embedding_obsm_key="X_pca",
    n_jobs=-1)

bm.benchmark()
    
bm.plot_results_table(min_max_scale=False, save_dir=os.path.join(PATH, "results/benchmark/combined/celltype"))
df = bm.get_results(min_max_scale=False)
df.to_csv(os.path.join(PATH, "results/benchmark/combined/celltype/bm_combined.csv"))

#
# bm = Benchmarker(
#     adata_sub,
#     batch_key="batch",
#     label_key="subtype_celltype",
#     bio_conservation_metrics=biocons,
#     batch_correction_metrics=batchcons,
#     embedding_obsm_keys=["Unintegrated", "X_scANVI_ct", "X_scANVI_sub_ct", "X_scVI", "X_seurat_harmony", "X_seurat_rpca", "X_seurat_mnn"],
#     pre_integrated_embedding_obsm_key="X_pca",
#     n_jobs=-1)
#
# bm.benchmark()
#
# bm.plot_results_table(min_max_scale=False, save_dir=os.path.join(PATH, "results/benchmark/combined/subtype_celltype"))
# df = bm.get_results(min_max_scale=False)
# df.to_csv(os.path.join(PATH, "results/benchmark/combined/subtype_celltype/bm_combined.csv"))
#
# bm = Benchmarker(
#     adata_sub,
#     batch_key="batch",
#     label_key="celltype_new",
#     bio_conservation_metrics=biocons,
#     batch_correction_metrics=batchcons,
#     embedding_obsm_keys=["Unintegrated", "X_scANVI_ct", "X_scANVI_sub_ct", "X_scVI", "X_seurat_harmony", "X_seurat_rpca", "X_seurat_mnn", "X_concat"],
#     pre_integrated_embedding_obsm_key="X_pca",
#     n_jobs=-1)
# 
# bm.benchmark()
# 
# bm.plot_results_table(min_max_scale=False, save_dir=os.path.join(PATH, "results/benchmark/combined/celltype_new"))
# df = bm.get_results(min_max_scale=False)
# df.to_csv(os.path.join(PATH, "results/benchmark/combined/celltype_new/bm_combined.csv"))

with open(os.path.join(PATH, "results/benchmark/combined/celltype/bm_object.pkl"), 'wb') as outp:
    pickle.dump(bm, outp, pickle.HIGHEST_PROTOCOL)
