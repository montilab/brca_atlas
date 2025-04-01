library(tidyverse)
library(Seurat)

PATH <- file.path(Sys.getenv("MLAB"), "projects/brcameta/brca_atlas/")

# # Adding slots to main object
combined_seurat <- readRDS(file.path(PATH, "data/sc/combined_seurat.rds"))
scvi <- read.csv(file.path(PATH, "data/embeddings/all/scvi30.csv"), header = FALSE)
rownames(scvi) <- colnames(combined_seurat)
colnames(scvi) <- paste0("scvi_", 1:30)
scanvi_ct <- read.csv(file.path(PATH, "data/embeddings/all/scanvi_30_celltypist.csv"), header = FALSE)
rownames(scanvi_ct) <- colnames(combined_seurat)
colnames(scanvi_ct) <- paste0("scanvi_ct_", 1:30)
scanvi_sub_ct <- read.csv(file.path(PATH, "data/embeddings/all/scanvi_30_subtype_celltypist.csv"), header = FALSE)
rownames(scanvi_sub_ct) <- colnames(combined_seurat)
colnames(scanvi_sub_ct) <- paste0("scanvi_sub_ct_", 1:30)
harmony <- read.csv(file.path(PATH, "data/embeddings/all/harmony_embedding_combined.csv"), header = TRUE, row.names = 1)
stopifnot(all.equal(rownames(harmony), colnames(combined_seurat)))
mnn <- read.csv(file.path(PATH, "data/embeddings/all/mnn_embedding_combined.csv"), header = TRUE, row.names = 1)
stopifnot(all.equal(rownames(mnn), colnames(combined_seurat)))
rpca <- read.csv(file.path(PATH, "data/embeddings/all/rpca_embedding_combined.csv"), header = TRUE, row.names = 1)
stopifnot(all.equal(rownames(rpca), colnames(combined_seurat)))

combined_seurat@reductions$harmony <- CreateDimReducObject(embeddings = as.matrix(harmony),
                                                               key = "harmony_",
                                                               assay = DefaultAssay(combined_seurat))
combined_seurat@reductions$integrated.mnn <- CreateDimReducObject(embeddings = as.matrix(mnn),
                                                                 key = "mnn_",
                                                                 assay = DefaultAssay(combined_seurat))
combined_seurat@reductions$integrated.rpca <- CreateDimReducObject(embeddings = as.matrix(rpca),
                                                            key = "integratedrpca_",
                                                            assay = DefaultAssay(combined_seurat))
combined_seurat@reductions$scvi <- CreateDimReducObject(embeddings = as.matrix(scvi),
                                                        key = "scvi_",
                                                        assay = DefaultAssay(combined_seurat))
combined_seurat@reductions$scanvi_ct <- CreateDimReducObject(embeddings = as.matrix(scanvi_ct),
                                                          key = "scanvi_ct_",
                                                          assay = DefaultAssay(combined_seurat))
combined_seurat@reductions$scanvi_sub_ct <- CreateDimReducObject(embeddings = as.matrix(scanvi_sub_ct),
                                                           key = "scanvi_sub_ct_",
                                                           assay = DefaultAssay(combined_seurat))
rm(harmony)
rm(mnn)
rm(rpca)
rm(scvi)
rm(scanvi_ct)
rm(scanvi_sub_ct)
gc()
saveRDS(combined_seurat, file.path(PATH, "data/sc/combined_seurat_filtered_int.rds"))

combined_seurat <- readRDS(file.path(PATH, "data/sc/combined_seurat_filtered_int.rds"))
# combined_seurat <- RunUMAP(combined_seurat, reduction = "pca", dims = 1:200, reduction.name = "umap.pca")
combined_seurat <- RunUMAP(combined_seurat, reduction = "scvi", dims = 1:30, reduction.name = "umap.scvi")
combined_seurat <- RunUMAP(combined_seurat, reduction = "scanvi_ct", dims = 1:30, reduction.name = "umap.scanvi_ct")
combined_seurat <- RunUMAP(combined_seurat, reduction = "scanvi_sub_ct", dims = 1:30, reduction.name = "umap.scanvi_sub_ct")
combined_seurat <- RunUMAP(combined_seurat, reduction = "harmony", dims = 1:200, reduction.name = "umap.harmony")
combined_seurat <- RunUMAP(combined_seurat, reduction = "integrated.rpca", dims = 1:200, reduction.name = "umap.rpca")
combined_seurat <- RunUMAP(combined_seurat, reduction = "integrated.mnn", dims = 1:50, reduction.name = "umap.mnn")
saveRDS(combined_seurat, file.path(PATH, "data/sc/combined_seurat_filtered_int.rds"))
