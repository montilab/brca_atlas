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
id <- as.numeric(Sys.getenv("SGE_TASK_ID"))
to_do <- c(17,23,24,25,26,27,28,29,30)
id <- to_do[id]
print(id)

Idents(diet_seurat_obj) <- diet_seurat_obj$RNA_snn_res.0.1
markers <- FindMarkers(diet_seurat_obj,
                       ident.1 = id,
                       assay = "RNA",
                       test.use = "MAST",
                       latent.vars = "batch",
                       only.pos = TRUE,
                       min.cells.group = 3,
                       verbose = TRUE)
saveRDS(markers, file.path(PATH,
                           paste0("results/clustering/",compartment,"_cluster_",id, "_01_mast_markers.rds")))
