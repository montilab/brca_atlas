library(tidyverse)
library(Seurat)

PATH <- file.path(Sys.getenv("MLAB"), "projects/brcameta/brca_atlas/")
do_save <- FALSE

# Seurat Obj Paths
paths <- list(epi = file.path(PATH, "data/sc/epi_rpca_subset.rds"),
              imm = file.path(PATH, "data/sc/imm_rpca_subset.rds"),
              strom = file.path(PATH, "data/sc/strom_rpca_subset.rds"))

# 1. Find All Markers (MAST)
compartment <- "epi"
diet_seurat_obj <- readRDS(paths[[compartment]])

Idents(diet_seurat_obj) <- diet_seurat_obj$RNA_snn_res.0.2
markers <- FindAllMarkers(diet_seurat_obj,
                          assay = "RNA",
                          test.use = "MAST",
                          latent.vars = "batch",
                          only.pos = TRUE,
                          min.cells.group = 3,
                          verbose = TRUE)
saveRDS(markers, file.path(PATH,
                           paste0("results/clustering/",compartment,"_0.2_mast_markers.rds")))
