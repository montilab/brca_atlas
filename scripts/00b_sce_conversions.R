library(Seurat)
library(tidyverse)
library(sceasy)
library(reticulate)
reticulate::use_condaenv("r-sceasy")
options(Seurat.object.assay.version = "v5")

PATH <- file.path(Sys.getenv("MLAB"), "projects/brcameta/brca_atlas/")
DATA_PATH <- file.path(Sys.getenv("CBM"), "otherStudies/scRNAseq/brca")

# Loading Data
wu_natgen_2021_sce <- readRDS(file=file.path(DATA_PATH, "wu_natgen_2021/processed/processed_sce.rds"))
wu_natgen_2021_sce$batch <- "wu_natgen_2021"
wu_embo_2020_sce <- readRDS(file=file.path(DATA_PATH, "wu_embo_2020/processed/processed_sce.rds"))
wu_embo_2020_sce$batch <- "wu_embo_2020"
wu_genomemed_2021_sce <- readRDS(file=file.path(DATA_PATH, "wu_genomemed_2021/processed_sce.rds"))
wu_genomemed_2021_sce$batch <- "wu_genomemed_2021"
pal_cancer <- readRDS(file=file.path(DATA_PATH, "pal_2021/processed/combined_cancer.rds"))
pal_cancer$batch <- "pal_2021"
pal_cancer$singler_pred <- pal_cancer$celltype_pred2
pal_cancer$celltype_pred2 <- NULL
qian_sce <- readRDS(file=file.path(DATA_PATH, "qian_2020/processed_sce_v2.rds"))
qian_sce$batch <- "qian_2020"
gao_sce <- readRDS(file=file.path(DATA_PATH, "gao_2021/processed_sce.rds"))
gao_sce$batch <- "gao_2021"
tietscher_combined <- readRDS(file=file.path(DATA_PATH, "tietscher_2023/processed/combined.rds"))
tietscher_combined$batch <- "tietscher_2023"
bassez_2021 <- readRDS(file=file.path(DATA_PATH, "bassez_2021/processed_cohort1.rds"))
bassez_2021$batch <- "bassez_2021"
xu_combined <- readRDS(file=file.path(DATA_PATH, "xu_2021/processed/combined.rds"))
xu_combined$batch <- "xu_2021"
barkley_combined <- readRDS(file=file.path(DATA_PATH, "barkley_2022/processed/combined.rds"))
barkley_combined$batch <- "barkley_2022"
liu_sce <- readRDS(file=file.path(DATA_PATH, "liu_2023/processed_sce.rds"))
liu_sce$batch <- "liu_2023"
wang_sce <- readRDS(file=file.path(DATA_PATH, "wang_2024/processed/ER/combined.rds"))
wang_sce$batch <- "wang_2024"
wang_sce2 <- readRDS(file=file.path(DATA_PATH, "wang_2024/processed/TNBC/combined.rds"))
wang_sce2$batch <- "wang_2024"
wang_combined <- list(wang_sce, wang_sce2)
test_rownames <- lapply(wang_combined, rownames)
test_intersect_rownames <- purrr::reduce(test_rownames, intersect)
wang_combined <- lapply(wang_combined, function(x) x[rownames(x) %in% test_intersect_rownames,])
wang_combined <- purrr::reduce(wang_combined, SingleCellExperiment::cbind)

# Making Cell IDs unique
make_unique <- function(x) {
  dups <- duplicated(x)
  if (any(dups)) {
    suffix <- seq_along(x)
    x[dups] <- paste0(x[dups], "_", suffix[dups])
  }
  return(x)
}
pal_ids <- make_unique(colnames(pal_cancer))
colnames(pal_cancer) <- pal_ids
rownames(pal_cancer@colData) <- pal_ids

xu_ids <- make_unique(colnames(xu_combined))
colnames(xu_combined) <- xu_ids
rownames(xu_combined@colData) <- xu_ids

barkley_ids <- make_unique(colnames(barkley_combined))
colnames(barkley_combined) <- barkley_ids
rownames(barkley_combined@colData) <- barkley_ids

wang_ids <- make_unique(colnames(wang_combined))
colnames(wang_combined) <- wang_ids
rownames(wang_combined@colData) <- wang_ids

# Converting to Seurat
wu_natgen_2021_sce <- as.Seurat(wu_natgen_2021_sce)
wu_embo_2020_sce <- as.Seurat(wu_embo_2020_sce)
wu_genomemed_2021_sce <- as.Seurat(wu_genomemed_2021_sce)
pal_cancer <- as.Seurat(pal_cancer)
qian_sce <- as.Seurat(qian_sce)
gao_sce <- as.Seurat(gao_sce)
tietscher_combined <- as.Seurat(tietscher_combined)
bassez_2021 <- as.Seurat(bassez_2021)
xu_combined <- as.Seurat(xu_combined)
barkley_combined <- as.Seurat(barkley_combined)
liu_sce <- as.Seurat(liu_sce)
wang_combined <- as.Seurat(wang_combined)

# Variable Feature Selection, Scaling, and PCA per dataset
wu_natgen_2021_sce <- FindVariableFeatures(wu_natgen_2021_sce, selection.method = "vst", nfeatures=3000)
wu_natgen_2021_sce <- ScaleData(wu_natgen_2021_sce)
wu_natgen_2021_sce <- RunPCA(wu_natgen_2021_sce, npcs = 50)
wu_natgen_2021_sce <- RunUMAP(wu_natgen_2021_sce,
                              reduction = "pca",
                              dims = 1:50,
                              reduction.name = "umap.pca")

wu_embo_2020_sce <- FindVariableFeatures(wu_embo_2020_sce, selection.method = "vst", nfeatures=3000)
wu_embo_2020_sce <- ScaleData(wu_embo_2020_sce)
wu_embo_2020_sce <- RunPCA(wu_embo_2020_sce, npcs = 50)
wu_embo_2020_sce <- RunUMAP(wu_embo_2020_sce,
                            reduction = "pca",
                            dims = 1:50,
                            reduction.name = "umap.pca")

wu_genomemed_2021_sce <- FindVariableFeatures(wu_genomemed_2021_sce, selection.method = "vst", nfeatures=3000)
wu_genomemed_2021_sce <- ScaleData(wu_genomemed_2021_sce)
wu_genomemed_2021_sce <- RunPCA(wu_genomemed_2021_sce, npcs = 50)
wu_genomemed_2021_sce <- RunUMAP(wu_genomemed_2021_sce,
                                 reduction = "pca",
                                 dims = 1:50,
                                 reduction.name = "umap.pca")

pal_cancer <- FindVariableFeatures(pal_cancer, selection.method = "vst", nfeatures=3000)
pal_cancer <- ScaleData(pal_cancer)
pal_cancer <- RunPCA(pal_cancer, npcs = 50)
pal_cancer <- RunUMAP(pal_cancer,
                      reduction = "pca",
                      dims = 1:50,
                      reduction.name = "umap.pca")

qian_sce <- FindVariableFeatures(qian_sce, selection.method = "vst", nfeatures=3000)
qian_sce <- ScaleData(qian_sce)
qian_sce <- RunPCA(qian_sce, npcs = 50)
qian_sce <- RunUMAP(qian_sce,
                    reduction = "pca",
                    dims = 1:50,
                    reduction.name = "umap.pca")

gao_sce <- FindVariableFeatures(gao_sce, selection.method = "vst", nfeatures=3000)
gao_sce <- ScaleData(gao_sce)
gao_sce <- RunPCA(gao_sce, npcs = 50)
gao_sce <- RunUMAP(gao_sce,
                   reduction = "pca",
                   dims = 1:50,
                   reduction.name = "umap.pca")

tietscher_combined <- FindVariableFeatures(tietscher_combined, selection.method = "vst", nfeatures=3000)
tietscher_combined <- ScaleData(tietscher_combined)
tietscher_combined <- RunPCA(tietscher_combined, npcs = 50)
tietscher_combined <- RunUMAP(tietscher_combined,
                              reduction = "pca",
                              dims = 1:50,
                              reduction.name = "umap.pca")

bassez_2021 <- FindVariableFeatures(bassez_2021, selection.method = "vst", nfeatures=3000)
bassez_2021 <- ScaleData(bassez_2021)
bassez_2021 <- RunPCA(bassez_2021, npcs = 50)
bassez_2021 <- RunUMAP(bassez_2021,
                       reduction = "pca",
                       dims = 1:50,
                       reduction.name = "umap.pca")

xu_combined <- FindVariableFeatures(xu_combined, selection.method = "vst", nfeatures=3000)
xu_combined <- ScaleData(xu_combined)
xu_combined <- RunPCA(xu_combined, npcs = 50)
xu_combined <- RunUMAP(xu_combined,
                       reduction = "pca",
                       dims = 1:50,
                       reduction.name = "umap.pca")

barkley_combined <- FindVariableFeatures(barkley_combined, selection.method = "vst", nfeatures=3000)
barkley_combined <- ScaleData(barkley_combined)
barkley_combined <- RunPCA(barkley_combined, npcs = 50)
barkley_combined <- RunUMAP(barkley_combined,
                            reduction = "pca",
                            dims = 1:50,
                            reduction.name = "umap.pca")

liu_sce <- FindVariableFeatures(liu_sce, selection.method = "vst", nfeatures=3000)
liu_sce <- ScaleData(liu_sce)
liu_sce <- RunPCA(liu_sce, npcs = 50)
liu_sce <- RunUMAP(liu_sce,
                   reduction = "pca",
                   dims = 1:50,
                   reduction.name = "umap.pca")

wang_combined <- FindVariableFeatures(wang_combined, selection.method = "vst", nfeatures=3000)
wang_combined <- ScaleData(wang_combined)
wang_combined <- RunPCA(wang_combined, npcs = 50)
wang_combined <- RunUMAP(wang_combined,
                         reduction = "pca",
                         dims = 1:50,
                         reduction.name = "umap.pca")

saveRDS(wu_natgen_2021_sce, file.path(DATA_PATH, "wu_natgen_2021/wu_natgen_seurat.rds"))
saveRDS(wu_embo_2020_sce, file.path(DATA_PATH, "wu_embo_2020/wu_embo_seurat.rds"))
saveRDS(wu_genomemed_2021_sce, file.path(DATA_PATH, "wu_genomemed_2021/wu_genomemed_seurat.rds"))
saveRDS(pal_cancer, file.path(DATA_PATH, "pal_2021/pal_combined_seurat.rds"))
saveRDS(qian_sce, file.path(DATA_PATH, "qian_2020/qian_seurat.rds"))
saveRDS(gao_sce, file.path(DATA_PATH, "gao_2021/gao_seurat.rds"))
saveRDS(tietscher_combined, file.path(DATA_PATH, "tietscher_2023/tietscher_combined_seurat.rds"))
saveRDS(bassez_2021, file.path(DATA_PATH, "bassez_2021/bassez_combined_seurat.rds"))
saveRDS(xu_combined, file.path(DATA_PATH, "xu_2021/xu_combined_seurat.rds"))
saveRDS(barkley_combined, file.path(DATA_PATH, "barkley_2022/barkley_combined_seurat.rds"))
saveRDS(liu_sce, file.path(DATA_PATH, "liu_2023/liu_seurat.rds"))
saveRDS(wang_combined, file.path(DATA_PATH, "wang_2024/wang_combined_seurat.rds"))

# Converting to AnnData
sceasy::convertFormat(wu_natgen_2021_sce, from="seurat", to="anndata",
                      outFile=file.path(DATA_PATH, "wu_natgen_2021/wu_natgen.h5ad"))
sceasy::convertFormat(wu_embo_2020_sce, from="seurat", to="anndata",
                      outFile=file.path(DATA_PATH, "wu_embo_2020/wu_embo.h5ad"))
sceasy::convertFormat(wu_genomemed_2021_sce, from="seurat", to="anndata",
                      outFile=file.path(DATA_PATH, "wu_genomemed_2021/wu_genomemed.h5ad"))
sceasy::convertFormat(pal_cancer, from="seurat", to="anndata",
                      outFile=file.path(DATA_PATH, "pal_2021/pal_combined.h5ad"))
sceasy::convertFormat(gao_sce, from="seurat", to="anndata",
                      outFile=file.path(DATA_PATH, "gao_2021/gao.h5ad"))
sceasy::convertFormat(tietscher_combined, from="seurat", to="anndata",
                      outFile=file.path(DATA_PATH, "tietscher_2023/tietscher_combined.h5ad"))
sceasy::convertFormat(bassez_2021, from="seurat", to="anndata",
                      outFile=file.path(DATA_PATH, "bassez_2021/bassez_combined.h5ad"))
sceasy::convertFormat(xu_combined, from="seurat", to="anndata",
                      outFile=file.path(DATA_PATH, "xu_2021/xu_combined.h5ad"))
sceasy::convertFormat(barkley_combined, from="seurat", to="anndata",
                      outFile=file.path(DATA_PATH, "barkley_2022/barkley_combined.h5ad"))
sceasy::convertFormat(liu_sce, from="seurat", to="anndata",
                      outFile=file.path(DATA_PATH, "liu_2023/liu.h5ad"))
sceasy::convertFormat(wang_combined, from="seurat", to="anndata",
                      outFile=file.path(DATA_PATH, "wang_2024/wang_combined.h5ad"))
sceasy::convertFormat(qian, from="seurat", to="anndata",
                      outFile=file.path(DATA_PATH, "qian_2020/qian.h5ad"))
