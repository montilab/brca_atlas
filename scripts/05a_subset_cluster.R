library(tidyverse)
library(Seurat)
library(scCustomize)
#library(reticulate)
#reticulate::use_condaenv("r-umap", required=TRUE)
#reticulate::py_config()
options(future.globals.maxSize = 100000 * 1024^2)
PATH <- file.path(Sys.getenv("MLAB"), "projects/brcameta/brca_atlas/")
source("util.R")
do_save <- TRUE

# combined_seurat_filtered <- readRDS(file.path(PATH, "data/sc/combined_seurat_filtered_int.rds"))
# 
# combined_seurat_filtered <- JoinLayers(combined_seurat_filtered)
# combined_seurat_filtered$celltypist_broad <- celltypist_broad(combined_seurat_filtered$celltypist_pred)

# 1. Comparing Annotations
# RPCA-RNA
combined_seurat_rpca <- readRDS(file.path(PATH, "data/sc/combined_seurat_rpca.rds"))
# combined_seurat_rpca <- DietSeurat(combined_seurat_filtered,
#                                    layers = c("counts", "data", "scale.data"),
#                                    assays = c("RNA"),
#                                    dimreducs = c("pca", "umap.pca",
#                                                  "integrated.rpca", "umap.rpca"))
# combined_seurat_rpca <- FindNeighbors(combined_seurat_rpca, dims = 1:200, reduction = "integrated.rpca")
# combined_seurat_rpca <- FindClusters(combined_seurat_rpca, resolution = 0.2)
# saveRDS(combined_seurat_rpca, file.path(PATH, "data/sc/combined_seurat_rpca.rds"))

# Comparing broad annotations
# p1 <- DimPlot_scCustom(combined_seurat_rpca,
#                        reduction = "umap.rpca",
#                        group.by = "RNA_snn_res.0.2",
#                        label.box = TRUE,
#                        repel = TRUE,
#                        raster = FALSE) + labs(title = "Clustering (0.2)", x = "UMAP 1", y = "UMAP 2") +
#   theme(legend.position="none")
# p2 <- DimPlot_scCustom(combined_seurat_rpca,
#                        reduction = "umap.rpca",
#                        group.by = "singleR_broad",
#                        raster = FALSE) + labs(title = "SingleR", x = "UMAP 1", y = "UMAP 2")
# p3 <- DimPlot_scCustom(combined_seurat_rpca,
#                        reduction = "umap.rpca",
#                        group.by = "celltypist_broad",
#                        raster = FALSE) + labs(title = "Celltypist", x = "UMAP 1", y = "UMAP 2")
# p4 <- DimPlot_scCustom(combined_seurat_rpca,
#                        reduction = "umap.rpca",
#                        group.by = "author_broad",
#                        raster = FALSE) + labs(title = "Author", x = "UMAP 1", y = "UMAP 2")
# combined <- cowplot::plot_grid(p1, p2, p3, p4, ncol=2, nrow = 2)
# 
# ggsave(plot = combined, filename = file.path(PATH, "results/umaps/all_clustering_broad_02.png"),
#        width = 9, height = 8)
# p1 <- DimPlot_scCustom(combined_seurat_rpca,
#                       reduction = "umap.rpca",
#                       group.by = "cluster_broad",
#                       label.box = TRUE,
#                       raster = FALSE)

combined_seurat_rpca$cluster_broad <- with(combined_seurat_rpca@meta.data,
                                           case_when(`RNA_snn_res.0.2` == 0 ~ "Immune",
                                                     `RNA_snn_res.0.2` == 1 ~ "Epithelial",
                                                     `RNA_snn_res.0.2` == 2 ~ "Immune",
                                                     `RNA_snn_res.0.2` == 3 ~ "Immune",
                                                     `RNA_snn_res.0.2` == 4 ~ "Stromal",
                                                     `RNA_snn_res.0.2` == 5 ~ "Epithelial",
                                                     `RNA_snn_res.0.2` == 6 ~ "Stromal",
                                                     `RNA_snn_res.0.2` == 7 ~ "Epithelial",
                                                     `RNA_snn_res.0.2` == 8 ~ "Immune",
                                                     `RNA_snn_res.0.2` == 9 ~ "Stromal",
                                                     `RNA_snn_res.0.2` == 10 ~ "Immune",
                                                     `RNA_snn_res.0.2` == 11 ~ "Immune",
                                                     `RNA_snn_res.0.2` == 12 ~ "Epithelial",
                                                     `RNA_snn_res.0.2` == 13 ~ "Epithelial",
                                                     `RNA_snn_res.0.2` == 14 ~ "Epithelial",
                                                     `RNA_snn_res.0.2` == 15 ~ "Epithelial",
                                                     `RNA_snn_res.0.2` == 16 ~ "Epithelial",
                                                     `RNA_snn_res.0.2` == 17 ~ "Epithelial",
                                                     `RNA_snn_res.0.2` == 18 ~ "Epithelial",
                                                     `RNA_snn_res.0.2` == 19 ~ "Epithelial",
                                                     `RNA_snn_res.0.2` == 20 ~ "Epithelial",
                                                     `RNA_snn_res.0.2` == 21 ~ "Epithelial",
                                                     `RNA_snn_res.0.2` == 22 ~ "Immune",
                                                     `RNA_snn_res.0.2` == 23 ~ "Epithelial",
                                                     `RNA_snn_res.0.2` == 24 ~ "Immune",
                                                     `RNA_snn_res.0.2` == 25 ~ "Epithelial",
                                                     `RNA_snn_res.0.2` == 26 ~ "Immune",
                                                     `RNA_snn_res.0.2` == 27 ~ "Epithelial",
                                                     `RNA_snn_res.0.2` == 28 ~ "Epithelial",
                                                     `RNA_snn_res.0.2` == 29 ~ "Epithelial",
                                                     `RNA_snn_res.0.2` == 30 ~ "Epithelial",
                                                     `RNA_snn_res.0.2` == 31 ~ "Stromal",
                                                     `RNA_snn_res.0.2` == 32 ~ "Stromal",
                                                     `RNA_snn_res.0.2` == 33 ~ "Epithelial"))


## Clustering each broad group
epi_subset <- subset(combined_seurat_rpca,
                     (cluster_broad == "Epithelial"))
epi_subset <- FindNeighbors(epi_subset, dims = 1:200, reduction = "integrated.rpca")
for (res in c(0.1, 0.2, 0.4, 0.6, 0.8)) {
  epi_subset <- FindClusters(epi_subset, resolution = res)
}
# imm_subset <- subset(combined_seurat_rpca,
#                      (cluster_broad == "Immune"))
# imm_subset <- FindNeighbors(imm_subset, dims = 1:200, reduction = "integrated.rpca")
# for (res in c(0.2, 0.4, 0.6, 0.8)) {
#   imm_subset <- FindClusters(imm_subset, resolution = res)
# }
strom_subset <- subset(combined_seurat_rpca,
                       (cluster_broad == "Stromal"))
strom_subset <- FindNeighbors(strom_subset, dims = 1:200, reduction = "integrated.rpca")
for (res in c(0.2, 0.4, 0.6, 0.8)) {
  strom_subset <- FindClusters(strom_subset, resolution = res)
}

# imm_subset$celltypist_shortened <- celltypist_short(imm_subset$celltypist_pred)
strom_subset$celltypist_shortened <- celltypist_short(strom_subset$celltypist_pred)
epi_subset$celltypist_shortened <- celltypist_short(epi_subset$celltypist_pred)


# imm_subset <- RunUMAP(imm_subset, reduction = "integrated.rpca", dims = 1:200, reduction.name = "umap.rpca")
epi_subset <- RunUMAP(epi_subset, reduction = "integrated.rpca", dims = 1:200, reduction.name = "umap.rpca")
strom_subset <- RunUMAP(strom_subset, reduction = "integrated.rpca", dims = 1:200, reduction.name = "umap.rpca")

saveRDS(epi_subset, file.path(PATH, "data/sc/epi_rpca_subset.rds"))
# saveRDS(imm_subset, file.path(PATH, "data/sc/imm_rpca_subset.rds"))
saveRDS(strom_subset, file.path(PATH, "data/sc/strom_rpca_subset.rds"))
saveRDS(combined_seurat_rpca, file.path(PATH, "data/sc/combined_seurat_rpca.rds"))
