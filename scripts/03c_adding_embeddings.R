library(tidyverse)
library(Seurat)

PATH <- file.path(Sys.getenv("MLAB"), "projects/brcameta/brca_atlas/")

# # Adding slots to main object
combined_seurat <- readRDS(file.path(PATH, "data/sc/combined_seurat_filtered_int.rds"))
scvi <- read.csv(file.path(PATH, "data/embeddings/all/scvi.csv"), header = FALSE)
rownames(scvi) <- colnames(combined_seurat)
colnames(scvi) <- paste0("scvi_", 1:200)
scanvi_ct <- read.csv(file.path(PATH, "data/embeddings/all/scanvi_celltypist.csv"), header = FALSE)
rownames(scanvi_ct) <- colnames(combined_seurat)
colnames(scanvi_ct) <- paste0("scanvi_ct_", 1:200)
scanvi_sub_ct <- read.csv(file.path(PATH, "data/embeddings/all/scanvi_subtype_celltypist.csv"), header = FALSE)
rownames(scanvi_sub_ct) <- colnames(combined_seurat)
colnames(scanvi_sub_ct) <- paste0("scanvi_sub_ct_", 1:200)
# harmony <- read.csv(file.path(PATH, "data/embeddings/all/harmony_embedding_combined.csv"), header = TRUE, row.names = 1)
# stopifnot(all.equal(rownames(harmony), colnames(combined_seurat)))
# # mnn <- read.csv(file.path(PATH, "data/embeddings/all/mnn_embedding_combined.csv"), header = TRUE, row.names = 1)
# # stopifnot(all.equal(rownames(mnn), colnames(combined_seurat)))
# # rpca <- read.csv(file.path(PATH, "data/embeddings/all/rpca_embedding_combined.csv"), header = TRUE, row.names = 1)
# # stopifnot(all.equal(rownames(rpca), colnames(combined_seurat)))
# 
# combined_seurat@reductions$harmony <- CreateDimReducObject(embeddings = as.matrix(harmony),
#                                                                key = "harmony_",
#                                                                assay = DefaultAssay(combined_seurat))
# # combined_seurat@reductions$integrated.mnn <- CreateDimReducObject(embeddings = as.matrix(mnn),
# #                                                                  key = "mnn_",
# #                                                                  assay = DefaultAssay(combined_seurat))
# # combined_seurat@reductions$integrated.rpca <- CreateDimReducObject(embeddings = as.matrix(rpca),
# #                                                             key = "integratedrpca_",
# #                                                             assay = DefaultAssay(combined_seurat))
combined_seurat@reductions$scvi <- CreateDimReducObject(embeddings = as.matrix(scvi),
                                                        key = "scvi_",
                                                        assay = DefaultAssay(combined_seurat))
combined_seurat@reductions$scanvi_ct <- CreateDimReducObject(embeddings = as.matrix(scanvi_ct),
                                                          key = "scanvi_ct_",
                                                          assay = DefaultAssay(combined_seurat))
combined_seurat@reductions$scanvi_sub_ct <- CreateDimReducObject(embeddings = as.matrix(scanvi_sub_ct),
                                                           key = "scanvi_sub_ct_",
                                                           assay = DefaultAssay(combined_seurat))
# rm(harmony)
# # rm(mnn)
# # rm(rpca)
rm(scvi)
rm(scanvi_ct)
rm(scanvi_sub_ct)
gc()
# saveRDS(combined_seurat, file.path(PATH, "data/sc/combined_seurat_filtered_int.rds"))
# 
# # Adding slots to epi data
# combined_epi_subset <- readRDS(file.path(PATH, "data/sc/combined_epi_subset_int.rds"))
# scanvi_ct <- read.csv(file.path(PATH, "data/embeddings/epi/scanvi_celltypist.csv"), header = FALSE)
# rownames(scanvi_ct) <- colnames(combined_epi_subset)
# colnames(scanvi_ct) <- paste0("scanvi_ct_", 1:200)
# scanvi_sub_ct <- read.csv(file.path(PATH, "data/embeddings/epi/scanvi_subtype_celltypist.csv"), header = FALSE)
# rownames(scanvi_sub_ct) <- colnames(combined_epi_subset)
# colnames(scanvi_sub_ct) <- paste0("scanvi_sub_ct_", 1:200)
# scvi <- read.csv(file.path(PATH, "data/embeddings/epi/scvi.csv"), header = FALSE)
# rownames(scvi) <- colnames(combined_epi_subset)
# colnames(scvi) <- paste0("scvi_", 1:200)
# harmony <- read.csv(file.path(PATH, "data/embeddings/epi/harmony_embedding_epi.csv"), header = TRUE, row.names = 1)
# stopifnot(all.equal(rownames(harmony), colnames(combined_epi_subset)))
# mnn <- read.csv(file.path(PATH, "data/embeddings/epi/mnn_embedding_epi.csv"), header = TRUE, row.names = 1)
# stopifnot(all.equal(rownames(mnn), colnames(combined_epi_subset)))
# rpca <- read.csv(file.path(PATH, "data/embeddings/epi/rpca_embedding_epi.csv"), header = TRUE, row.names = 1)
# stopifnot(all.equal(rownames(rpca), colnames(combined_epi_subset)))
# 
# combined_epi_subset@reductions$harmony <- CreateDimReducObject(embeddings = as.matrix(harmony),
#                                                                   key = "harmony_",
#                                                                   assay = DefaultAssay(combined_epi_subset))
# combined_epi_subset@reductions$integrated.mnn <- CreateDimReducObject(embeddings = as.matrix(mnn),
#                                                                  key = "mnn_",
#                                                                  assay = DefaultAssay(combined_epi_subset))
# combined_epi_subset@reductions$integrated.rpca <- CreateDimReducObject(embeddings = as.matrix(rpca),
#                                                                key = "integratedrpca_",
#                                                                assay = DefaultAssay(combined_epi_subset))
# combined_epi_subset@reductions$scanvi_ct <- CreateDimReducObject(embeddings = as.matrix(scanvi_ct),
#                                                                  key = "scanvi_ct_",
#                                                                  assay = DefaultAssay(combined_epi_subset))
# combined_epi_subset@reductions$scanvi_sub_ct <- CreateDimReducObject(embeddings = as.matrix(scanvi_sub_ct),
#                                                                      key = "scanvi_sub_ct_",
#                                                                      assay = DefaultAssay(combined_epi_subset))
# combined_epi_subset@reductions$scvi <- CreateDimReducObject(embeddings = as.matrix(scvi),
#                                                             key = "scvi_",
#                                                             assay = DefaultAssay(combined_epi_subset))
# rm(harmony)
# rm(mnn)
# rm(rpca)
# rm(scvi)
# rm(scanvi_ct)
# rm(scanvi_sub_ct)
# gc()
# saveRDS(combined_epi_subset, file.path(PATH, "data/sc/combined_epi_subset_int.rds"))
# 
# # Adding slots to immune data
# combined_immune_subset <- readRDS(file.path(PATH, "data/sc/combined_immune_subset_int.rds"))
# scanvi_ct <- read.csv(file.path(PATH, "data/embeddings/imm/scanvi_celltypist.csv"), header = FALSE)
# rownames(scanvi_ct) <- colnames(combined_immune_subset)
# colnames(scanvi_ct) <- paste0("scanvi_ct_", 1:200)
# scanvi_sub_ct <- read.csv(file.path(PATH, "data/embeddings/imm/scanvi_subtype_celltypist.csv"), header = FALSE)
# rownames(scanvi_sub_ct) <- colnames(combined_immune_subset)
# colnames(scanvi_sub_ct) <- paste0("scanvi_sub_ct_", 1:200)
# scvi <- read.csv(file.path(PATH, "data/embeddings/imm/scvi.csv"), header = FALSE)
# rownames(scvi) <- colnames(combined_immune_subset)
# colnames(scvi) <- paste0("scvi_", 1:200)
# harmony <- read.csv(file.path(PATH, "data/embeddings/imm/harmony_embedding_immune.csv"), header = TRUE, row.names = 1)
# stopifnot(all.equal(rownames(harmony), colnames(combined_immune_subset)))
# # mnn <- read.csv(file.path(PATH, "data/embeddings/imm/mnn_embedding_immune.csv"), header = TRUE, row.names = 1)
# # stopifnot(all.equal(rownames(mnn), colnames(combined_immune_subset)))
# rpca <- read.csv(file.path(PATH, "data/embeddings/imm/rpca_embedding_immune.csv"), header = TRUE, row.names = 1)
# stopifnot(all.equal(rownames(rpca), colnames(combined_immune_subset)))
# 
# combined_immune_subset@reductions$harmony <- CreateDimReducObject(embeddings = as.matrix(harmony),
#                                                                  key = "harmony_",
#                                                                  assay = DefaultAssay(combined_immune_subset))
# # combined_immune_subset@reductions$integrated.mnn <- CreateDimReducObject(embeddings = as.matrix(mnn),
# #                                                                  key = "mnn_",
# #                                                                  assay = DefaultAssay(combined_immune_subset))
# combined_immune_subset@reductions$integrated.rpca <- CreateDimReducObject(embeddings = as.matrix(rpca),
#                                                               key = "integratedrpca_",
#                                                               assay = DefaultAssay(combined_immune_subset))
# combined_immune_subset@reductions$scanvi_ct <- CreateDimReducObject(embeddings = as.matrix(scanvi_ct),
#                                                                     key = "scanvi_ct_",
#                                                                     assay = DefaultAssay(combined_immune_subset))
# combined_immune_subset@reductions$scanvi_sub_ct <- CreateDimReducObject(embeddings = as.matrix(scanvi_sub_ct),
#                                                                         key = "scanvi_sub_ct_",
#                                                                         assay = DefaultAssay(combined_immune_subset))
# combined_immune_subset@reductions$scvi <- CreateDimReducObject(embeddings = as.matrix(scvi),
#                                                                key = "scvi_",
#                                                                assay = DefaultAssay(combined_immune_subset))
# rm(harmony)
# # rm(mnn)
# rm(rpca)
# rm(scvi)
# rm(scanvi_ct)
# rm(scanvi_sub_ct)
# gc()
# saveRDS(combined_immune_subset, file.path(PATH, "data/sc/combined_immune_subset_int.rds"))
# 
# Adding slots to strom data
# combined_strom_subset <- readRDS(file.path(PATH, "data/sc/combined_strom_subset_int.rds"))
# scanvi_ct <- read.csv(file.path(PATH, "data/embeddings/strom/scanvi_celltypist.csv"), header = FALSE)
# rownames(scanvi_ct) <- colnames(combined_strom_subset)
# colnames(scanvi_ct) <- paste0("scanvi_ct_", 1:200)
# scanvi_sub_ct <- read.csv(file.path(PATH, "data/embeddings/strom/scanvi_subtype_celltypist.csv"), header = FALSE)
# rownames(scanvi_sub_ct) <- colnames(combined_strom_subset)
# colnames(scanvi_sub_ct) <- paste0("scanvi_sub_ct_", 1:200)
# scvi <- read.csv(file.path(PATH, "data/embeddings/strom/scvi.csv"), header = FALSE)
# rownames(scvi) <- colnames(combined_strom_subset)
# colnames(scvi) <- paste0("scvi_", 1:200)
# harmony <- read.csv(file.path(PATH, "data/embeddings/strom/harmony_embedding_strom.csv"), header = TRUE, row.names = 1)
# stopifnot(all.equal(rownames(harmony), colnames(combined_strom_subset)))
# mnn <- read.csv(file.path(PATH, "data/embeddings/strom/mnn_embedding_strom.csv"), header = TRUE, row.names = 1)
# stopifnot(all.equal(rownames(mnn), colnames(combined_strom_subset)))
# rpca <- read.csv(file.path(PATH, "data/embeddings/strom/rpca_embedding_strom.csv"), header = TRUE, row.names = 1)
# stopifnot(all.equal(rownames(rpca), colnames(combined_strom_subset)))
# 
# combined_strom_subset@reductions$harmony <- CreateDimReducObject(embeddings = as.matrix(harmony),
#                                                                  key = "harmony_",
#                                                                  assay = DefaultAssay(combined_strom_subset))
# combined_strom_subset@reductions$integrated.mnn <- CreateDimReducObject(embeddings = as.matrix(mnn),
#                                                                  key = "mnn_",
#                                                                  assay = DefaultAssay(combined_strom_subset))
# combined_strom_subset@reductions$integrated.rpca <- CreateDimReducObject(embeddings = as.matrix(rpca),
#                                                                  key = "integratedrpca_",
#                                                                  assay = DefaultAssay(combined_strom_subset))
# 
# combined_strom_subset@reductions$scanvi_ct <- CreateDimReducObject(embeddings = as.matrix(scanvi_ct),
#                                                                  key = "scanvi_ct_",
#                                                                  assay = DefaultAssay(combined_strom_subset))
# combined_strom_subset@reductions$scanvi_sub_ct <- CreateDimReducObject(embeddings = as.matrix(scanvi_sub_ct),
#                                                                      key = "scanvi_sub_ct_",
#                                                                      assay = DefaultAssay(combined_strom_subset))
# combined_strom_subset@reductions$scvi <- CreateDimReducObject(embeddings = as.matrix(scvi),
#                                                             key = "scvi_",
#                                                             assay = DefaultAssay(combined_strom_subset))
# rm(harmony)
# rm(mnn)
# rm(rpca)
# rm(scvi)
# rm(scanvi_ct)
# rm(scanvi_sub_ct)
# gc()
# saveRDS(combined_strom_subset, file.path(PATH, "data/sc/combined_strom_subset_int.rds"))

# Running UMAPs

#combined_strom_subset <- readRDS(file.path(PATH, "data/sc/combined_strom_subset_int.rds"))
# combined_strom_subset <- RunUMAP(combined_strom_subset, reduction = "scvi", dims = 1:200, reduction.name = "umap.scvi")
# combined_strom_subset <- RunUMAP(combined_strom_subset, reduction = "scanvi_sub_ct", dims = 1:200, reduction.name = "umap.scanvi_sub_ct")
# combined_strom_subset <- RunUMAP(combined_strom_subset, reduction = "scanvi_ct", dims = 1:200, reduction.name = "umap.scanvi_ct")
# combined_strom_subset <- RunUMAP(combined_strom_subset, reduction = "harmony", dims = 1:200, reduction.name = "umap.harmony")
# combined_strom_subset <- RunUMAP(combined_strom_subset, reduction = "integrated.rpca", dims = 1:200, reduction.name = "umap.rpca")
#combined_strom_subset <- RunUMAP(combined_strom_subset, reduction = "integrated.mnn", dims = 1:50, reduction.name = "umap.mnn")
# saveRDS(combined_strom_subset, file.path(PATH, "data/sc/combined_strom_subset_int.rds"))
# rm(combined_strom_subset)
# gc()
# 
# combined_epi_subset <- readRDS(file.path(PATH, "data/sc/combined_epi_subset_int.rds"))
# combined_epi_subset <- RunUMAP(combined_epi_subset, reduction = "scvi", dims = 1:200, reduction.name = "umap.scvi")
# combined_epi_subset <- RunUMAP(combined_epi_subset, reduction = "scanvi_sub_ct", dims = 1:200, reduction.name = "umap.scanvi_sub_ct")
# combined_epi_subset <- RunUMAP(combined_epi_subset, reduction = "scanvi_ct", dims = 1:200, reduction.name = "umap.scanvi_ct")
# combined_epi_subset <- RunUMAP(combined_epi_subset, reduction = "harmony", dims = 1:200, reduction.name = "umap.harmony")
# combined_epi_subset <- RunUMAP(combined_epi_subset, reduction = "integrated.rpca", dims = 1:200, reduction.name = "umap.rpca")
# # combined_epi_subset <- RunUMAP(combined_epi_subset, reduction = "integrated.mnn", dims = 1:50, reduction.name = "umap.mnn")
# saveRDS(combined_epi_subset, file.path(PATH, "data/sc/combined_epi_subset_int.rds"))
# rm(combined_epi_subset)
# gc()
# 
# # combined_immune_subset <- readRDS(file.path(PATH, "data/sc/combined_immune_subset_int.rds"))
# combined_immune_subset <- RunUMAP(combined_immune_subset, reduction = "scvi", dims = 1:200, reduction.name = "umap.scvi")
# combined_immune_subset <- RunUMAP(combined_immune_subset, reduction = "scanvi_sub_ct", dims = 1:200, reduction.name = "umap.scanvi_sub_ct")
# combined_immune_subset <- RunUMAP(combined_immune_subset, reduction = "scanvi_ct", dims = 1:200, reduction.name = "umap.scanvi_ct")
# combined_immune_subset <- RunUMAP(combined_immune_subset, reduction = "harmony", dims = 1:200, reduction.name = "umap.harmony")
# combined_immune_subset <- RunUMAP(combined_immune_subset, reduction = "integrated.rpca", dims = 1:200, reduction.name = "umap.rpca")
# # combined_immune_subset <- RunUMAP(combined_immune_subset, reduction = "integrated.mnn", dims = 1:50, reduction.name = "umap.mnn")
# saveRDS(combined_immune_subset, file.path(PATH, "data/sc/combined_immune_subset_int.rds"))
# rm(combined_immune_subset)
# gc()

# combined_seurat <- readRDS(file.path(PATH, "data/sc/combined_seurat_filtered_int.rds"))
# combined_seurat <- RunUMAP(combined_seurat, reduction = "pca", dims = 1:200, reduction.name = "umap.pca")
combined_seurat <- RunUMAP(combined_seurat, reduction = "scvi", dims = 1:200, reduction.name = "umap.scvi")
combined_seurat <- RunUMAP(combined_seurat, reduction = "scanvi_ct", dims = 1:200, reduction.name = "umap.scanvi_ct")
combined_seurat <- RunUMAP(combined_seurat, reduction = "scanvi_sub_ct", dims = 1:200, reduction.name = "umap.scanvi_sub_ct")
# combined_seurat <- RunUMAP(combined_seurat, reduction = "harmony", dims = 1:200, reduction.name = "umap.harmony")
# # combined_seurat <- RunUMAP(combined_seurat, reduction = "integrated.rpca", dims = 1:200, reduction.name = "umap.rpca")
# # combined_seurat <- RunUMAP(combined_seurat, reduction = "integrated.mnn", dims = 1:50, reduction.name = "umap.mnn")
saveRDS(combined_seurat, file.path(PATH, "data/sc/combined_seurat_filtered_int.rds"))

# imm_seurat <- readRDS(file.path(PATH, "data/sc/imm_seurat_diet.rds"))
# epi_seurat <- readRDS(file.path(PATH, "data/sc/epi_seurat_diet.rds"))
# strom_seurat <- readRDS(file.path(PATH, "data/sc/strom_seurat_diet.rds"))
# epi_seurat <- RunUMAP(epi_seurat, reduction = "scanvi_ct", dims = 1:200,  reduction.name = "umap.scanvi_ct")
# imm_seurat <- RunUMAP(imm_seurat, reduction = "scanvi_ct", dims = 1:200, reduction.name = "umap.scanvi_ct")
# strom_seurat <- RunUMAP(strom_seurat, reduction = "scanvi_ct", dims = 1:200, reduction.name = "umap.scanvi_ct")
# saveRDS(imm_seurat, file.path(PATH, "data/sc/imm_seurat_diet.rds"))
# saveRDS(epi_seurat, file.path(PATH, "data/sc/epi_seurat_diet.rds"))
# saveRDS(strom_seurat, file.path(PATH, "data/sc/strom_seurat_diet.rds"))