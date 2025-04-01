library(Seurat)
library(SeuratData)
library(SeuratWrappers)
library(doFuture)
options(future.globals.maxSize = Inf)

PATH <- file.path(Sys.getenv("MLAB"), "projects/brcameta/brca_atlas/")
seurat_obj <- readRDS(file.path(PATH, "data/sc/combined_seurat.rds"))
seurat_obj[["RNA"]] <- split(seurat_obj[["RNA"]], f = seurat_obj$batch)
# A consensus set is found after running find variable features on each dataset
seurat_obj <- FindVariableFeatures(seurat_obj, nfeatures=5000)
write.csv(VariableFeatures(seurat_obj), file.path(PATH, "data/var_genes/atlas_5000.csv"))
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj, npcs = 200)
seurat_obj <- RunUMAP(seurat_obj, reduction = "pca", dims = 1:200, reduction.name = "umap.pca")
saveRDS(seurat_obj, file.path(PATH, "data/sc/combined_seurat.rds"))
print(seurat_obj)


# 1. Harmony
seurat_obj <- IntegrateLayers(
  object = seurat_obj, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = TRUE)
print("Done with Harmony")
write.csv(seurat_obj@reductions$harmony@cell.embeddings,
          file.path(PATH, "data/embeddings/all/harmony_embedding_combined.csv"))

# 2. fastMNN
DefaultAssay(seurat_obj) <- "RNA"
seurat_obj <- IntegrateLayers(
  object = seurat_obj, method = FastMNNIntegration,
  new.reduction = "integrated.mnn",
  verbose = TRUE)
print("Done with fastMNN")
write.csv(seurat_obj@reductions$integrated.mnn@cell.embeddings,
          file.path(PATH, "data/embeddings/all/mnn_embedding_combined.csv"))

# 3. RPCA
# Use sctransform to get batch corrected counts through v5 RPCA. Otherwise only embedding is returned.
# seurat_obj <- SCTransform(seurat_obj)
# seurat_obj <- RunPCA(seurat_obj, npcs = 200)
# print(seurat_obj)
# seurat_obj <- IntegrateLayers(
#   object = seurat_obj, method = RPCAIntegration,
#   orig.reduction = "pca", new.reduction = "integrated.rpca", normalization.method = "SCT", k.anchor = 20,
#   verbose = TRUE)
seurat_obj <- IntegrateLayers(
  object = seurat_obj, method = RPCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.rpca", k.anchor = 20,
  verbose = TRUE)
print("Done with RPCA for")
write.csv(seurat_obj@reductions$integrated.rpca@cell.embeddings,
          file.path(PATH, paste0("data/embeddings/all/rpca_embedding_combined.csv")))

seurat_obj <- RunUMAP(seurat_obj, reduction = "harmony", dims = 1:200, reduction.name = "umap.harmony")
seurat_obj <- RunUMAP(seurat_obj, reduction = "integrated.mnn", dims = 1:50, reduction.name = "umap.mnn")
seurat_obj <- RunUMAP(seurat_obj, reduction = "integrated.rpca", dims = 1:200, reduction.name = "umap.rpca")

saveRDS(seurat_obj, file.path(PATH, "data/sc/combined_seurat_int.rds"))
