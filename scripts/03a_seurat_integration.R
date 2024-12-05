library(Seurat)
library(SeuratData)
library(SeuratWrappers)
library(doFuture)
# doFuture::registerDoFuture()
# future::plan("multisession", workers = 4)
# 
# # set this option when analyzing large datasets
# options(future.globals.maxSize = 30000 * 1024^2)
# memory.limit(size=20000)

PATH <- file.path(Sys.getenv("MLAB"), "projects/brcameta/brca_atlas/")

seurat_paths <- list(combined = list(read = file.path(PATH, "data/sc/combined_seurat_filtered.rds"),
                                     save = file.path(PATH, "data/sc/combined_seurat_filtered_int.rds")),
                     immune = list(read = file.path(PATH, "data/sc/combined_immune_subset.rds"),
                                   save = file.path(PATH, "data/sc/combined_immune_subset_int.rds")),
                     epi = list(read = file.path(PATH, "data/sc/combined_epi_subset.rds"),
                                save = file.path(PATH, "data/sc/combined_epi_subset_int.rds")),
                     strom = list(read = file.path(PATH, "data/sc/combined_strom_subset.rds"),
                                  save = file.path(PATH, "data/sc/combined_strom_subset_int.rds"))
                     )


foreach(seurat_path_name = names(seurat_paths)) %do% {
  seurat_path <- seurat_paths[[seurat_path_name]]
  seurat_obj <- readRDS(seurat_path[["read"]])
  print(seurat_obj)
  # 1. Harmony
  seurat_obj <- IntegrateLayers(
    object = seurat_obj, method = HarmonyIntegration,
    orig.reduction = "pca", new.reduction = "harmony",
    verbose = TRUE)
  print(paste0("Done with Harmony for ", seurat_path_name))
  seurat_obj <- RunUMAP(seurat_obj, reduction = "harmony", dims = 1:200, reduction.name = "umap.harmony")
  saveRDS(seurat_obj, seurat_path[["save"]])
  write.csv(seurat_obj@reductions$harmony@cell.embeddings,
            file.path(PATH, paste0("data/embeddings/harmony_embedding_", seurat_path_name, ".csv")))

  # 2. fastMNN
  DefaultAssay(seurat_obj) <- "RNA"
  seurat_obj <- IntegrateLayers(
    object = seurat_obj, method = FastMNNIntegration,
    new.reduction = "integrated.mnn",
    verbose = TRUE)
  print(paste0("Done with fastMNN for ", seurat_path_name))
  seurat_obj <- RunUMAP(seurat_obj, reduction = "integrated.mnn", dims = 1:50, reduction.name = "umap.mnn")
  saveRDS(seurat_obj, seurat_path[["save"]])
  write.csv(seurat_obj@reductions$integrated.mnn@cell.embeddings,
            file.path(PATH, paste0("data/embeddings/mnn_embedding_", seurat_path_name, ".csv")))

  # 3. RPCA
  # Using sctransform to get batch corrected counts through v5 RPCA
  seurat_obj <- SCTransform(seurat_obj)
  seurat_obj <- RunPCA(seurat_obj, npcs = 200)
  print(seurat_obj)
  seurat_obj <- IntegrateLayers(
    object = seurat_obj, method = RPCAIntegration,
    orig.reduction = "pca", new.reduction = "integrated.rpca", normalization.method = "SCT", k.anchor = 20,
    verbose = TRUE)
  print(paste0("Done with RPCA for ", seurat_path_name))
  seurat_obj <- RunUMAP(seurat_obj, reduction = "integrated.rpca", dims = 1:200, reduction.name = "umap.rpca")
  saveRDS(seurat_obj, seurat_path[["save"]])
  write.csv(seurat_obj@reductions$integrated.rpca@cell.embeddings,
            file.path(PATH, paste0("data/embeddings/rpca_embedding_", seurat_path_name, ".csv")))

}

