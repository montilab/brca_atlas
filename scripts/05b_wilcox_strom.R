library(tidyverse)
library(Seurat)

PATH <- file.path(Sys.getenv("MLAB"), "projects/brcameta/brca_atlas/")
do_save <- FALSE

# Seurat Obj Paths
paths <- list(epi = file.path(PATH, "data/sc/epi_rpca_subset.rds"),
              imm = file.path(PATH, "data/sc/imm_rpca_subset.rds"),
              strom = file.path(PATH, "data/sc/strom_rpca_subset.rds"))

# Find All Markers (Wilcox)
compartment <- "strom"
diet_seurat_obj <- readRDS(paths[[compartment]])

Idents(diet_seurat_obj) <- diet_seurat_obj$RNA_snn_res.0.4
markers <- FindAllMarkers(diet_seurat_obj,
                          assay = "RNA",
                          test.use = "wilcox",
                          only.pos = TRUE,
                          min.cells.group = 3,
                          verbose = TRUE)
saveRDS(markers, file.path(PATH,
                           paste0("results/clustering/",compartment,"_04_wilcox_markers.rds")))

# Find All Markers (wilcox with gene centering)
# 
diet_seurat_obj <- ScaleData(diet_seurat_obj, 
                             do.scale = FALSE, 
                             split.by = "batch",
                             features = rownames(diet_seurat_obj))
markers <- FindAllMarkers(diet_seurat_obj,
                          assay = "RNA",
                          test.use = "wilcox",
                          slot = "scale.data",
                          only.pos = TRUE,
                          min.cells.group = 3,
                          verbose = TRUE)
saveRDS(markers, file.path(PATH,
                           paste0("results/clustering/",compartment,"_04_wilcox_center_markers.rds")))

diet_seurat_obj <- ScaleData(diet_seurat_obj, 
                             split.by = "batch",
                             features = rownames(diet_seurat_obj))
markers <- FindAllMarkers(diet_seurat_obj,
                          assay = "RNA",
                          test.use = "wilcox",
                          slot = "scale.data",
                          only.pos = TRUE,
                          min.cells.group = 3,
                          verbose = TRUE)
saveRDS(markers, file.path(PATH,
                           paste0("results/clustering/",compartment,"_04_wilcox_scale_center_markers.rds")))
