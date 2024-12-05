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
# 
# 1. Comparing Annotations
# # Harmony
# combined_seurat_harmony <- DietSeurat(combined_seurat_filtered,
#                                       layers = c("counts", "data", "scale.data"),
#                                       assays = c("RNA"),
#                                       dimreducs = c("pca", "umap.pca",
#                                                     "harmony", "umap.harmony"))
# combined_seurat_harmony <- FindNeighbors(combined_seurat_harmony, reduction = "harmony")
# combined_seurat_harmony <- FindClusters(combined_seurat_harmony, resolution = 0.2)
# 
# saveRDS(combined_seurat_harmony, file.path(PATH, "data/sc/combined_seurat_harmony.rds"))
# # RPCA-RNA
# combined_seurat_rpca <- DietSeurat(combined_seurat_filtered,
#                                    layers = c("counts", "data", "scale.data"),
#                                    assays = c("RNA"),
#                                    dimreducs = c("pca", "umap.pca",
#                                                  "integrated.rpca", "umap.rpca"))
combined_seurat_rpca <- readRDS(file.path(PATH, "data/sc/combined_seurat_rpca.rds"))
# combined_seurat_rpca <- FindNeighbors(combined_seurat_rpca, dims = 1:200, reduction = "integrated.rpca")
#combined_seurat_rpca <- FindClusters(combined_seurat_rpca, method="igraph", graph.name = "RNA_snn", resolution = 0.1)
# 
# # RPCA-SCT
# DefaultAssay(combined_seurat_filtered) <- "SCT"
# combined_seurat_rpca <- DietSeurat(combined_seurat_filtered,
#                                    layers = c("counts", "data", "scale.data"),
#                                    assays = c("SCT"),
#                                    dimreducs = c("pca", "umap.pca",
#                                                  "integrated.rpca"))
# combined_seurat_rpca <- RunPCA(object = combined_seurat_rpca, 
#                                reduction.name = "pca2", 
#                                npcs = 200)
# combined_seurat_rpca <- RunUMAP(object = combined_seurat_rpca, reduction.name = "pca2")
# # Checking if pca of corrected counts is the same as rpca embedding
# comparison <- all.equal(combined_seurat_rpca@reductions$integrated.rpca@cell.embeddings, 
#                         combined_seurat_rpca@reductions$pca2@cell.embeddings)
# print(comparison)
# combined_seurat_rpca@reductions$pca2 <- NULL
# combined_seurat_rpca <- FindNeighbors(combined_seurat_rpca, reduction = "integrated.rpca")
# combined_seurat_rpca <- FindClusters(combined_seurat_rpca, graph.name = "RNA_snn", resolution = 0.2)
# 
# saveRDS(combined_seurat_rpca, file.path(PATH, "data/sc/combined_seurat_sct_rpca.rds"))
# 
# combined_seurat_scanvi_ct <- DietSeurat(combined_seurat_filtered,
#                                       layers = c("counts", "data", "scale.data"),
#                                       assays = c("RNA"),
#                                       dimreducs = c("pca", "umap.pca",
#                                                     "scanvi_ct", "umap.scanvi_ct"))
# combined_seurat_scanvi_ct <- FindNeighbors(combined_seurat_scanvi_ct, reduction = "scanvi_ct")
# combined_seurat_scanvi_ct <- FindClusters(combined_seurat_scanvi_ct, resolution = 0.2)
# 
# saveRDS(combined_seurat_scanvi_ct, file.path(PATH, "data/sc/combined_seurat_scanvi_ct.rds"))
# 
# p1 <- DimPlot_scCustom(combined_seurat_scanvi_ct, 
#                        reduction = "umap.scanvi_ct",
#                        group.by = "RNA_snn_res.0.2",
#                        label = TRUE,
#                        label.size = 2,
#                        label.box = TRUE, 
#                        repel = TRUE,
#                        raster = FALSE) + 
#                       labs(title = "Clustering Scanvi CT (0.2)", x = "UMAP 1", y = "UMAP 2") +
#                       theme(legend.position = "none")
# p2 <- DimPlot_scCustom(combined_seurat_scanvi_ct, 
#                        reduction = "umap.scanvi_ct",
#                        group.by = "singleR_broad",
#                        raster = FALSE) + labs(title = "SingleR", x = "UMAP 1", y = "UMAP 2")
# p3 <- DimPlot_scCustom(combined_seurat_scanvi_ct, 
#                        reduction = "umap.scanvi_ct",
#                        group.by = "celltypist_broad",
#                        raster = FALSE) + labs(title = "Celltypist", x = "UMAP 1", y = "UMAP 2")
# p4 <- DimPlot_scCustom(combined_seurat_rpca, 
#                        reduction = "umap.rpca",
#                        group.by = "RNA_snn_res.0.2",
#                        label = TRUE,
#                        label.box = TRUE, 
#                        repel = TRUE,
#                        raster = FALSE) + 
#                       labs(title = "Clustering RPCA (0.2)", x = "UMAP 1", y = "UMAP 2") +
#                       theme(legend.position = "none")
# p5 <- DimPlot_scCustom(combined_seurat_rpca, 
#                        reduction = "umap.rpca",
#                        group.by = "singleR_broad",
#                        raster = FALSE) + labs(title = "SingleR", x = "UMAP 1", y = "UMAP 2")
# p6 <- DimPlot_scCustom(combined_seurat_rpca, 
#                        reduction = "umap.rpca",
#                        group.by = "celltypist_broad",
#                        raster = FALSE) + labs(title = "Celltypist", x = "UMAP 1", y = "UMAP 2")
# p7 <- DimPlot_scCustom(combined_seurat_harmony, 
#                        reduction = "umap.harmony",
#                        group.by = "RNA_snn_res.0.2",
#                        label = TRUE,
#                        label.box = TRUE, 
#                        repel = TRUE,
#                        raster = FALSE) + 
#                       labs(title = "Clustering Harmony (0.2)", x = "UMAP 1", y = "UMAP 2") +
#                       theme(legend.position = "none")
# p8 <- DimPlot_scCustom(combined_seurat_harmony, 
#                        reduction = "umap.harmony",
#                        group.by = "singleR_broad",
#                        raster = FALSE) + labs(title = "SingleR", x = "UMAP 1", y = "UMAP 2")
# p9 <- DimPlot_scCustom(combined_seurat_harmony, 
#                        reduction = "umap.harmony",
#                        group.by = "celltypist_broad",
#                        raster = FALSE) + labs(title = "Celltypist", x = "UMAP 1", y = "UMAP 2")
# combined <- cowplot::plot_grid(p1, p2, p3,
#                                p4, p5, p6,
#                                p7, p8, p9, ncol = 3, nrow = 3)
# 
# ggsave(plot = combined, filename = file.path(PATH, "results/umaps/harmony_scanvi_rpca_comparison.png"),
#        width = 12, height = 10)
# 
# combined_seurat_harmony$cluster_broad <- with(combined_seurat_harmony@meta.data,
#                                            case_when(`RNA_snn_res.0.2` == 0 ~ "Immune",
#                                                      `RNA_snn_res.0.2` == 1 ~ "Epithelial",
#                                                      `RNA_snn_res.0.2` == 2 ~ "Immune",
#                                                      `RNA_snn_res.0.2` == 3 ~ "Epithelial",
#                                                      `RNA_snn_res.0.2` == 4 ~ "Epithelial",
#                                                      `RNA_snn_res.0.2` == 5 ~ "Stromal",
#                                                      `RNA_snn_res.0.2` == 6 ~ "Immune",
#                                                      `RNA_snn_res.0.2` == 7 ~ "Stromal",
#                                                      `RNA_snn_res.0.2` == 8 ~ "Immune",
#                                                      `RNA_snn_res.0.2` == 9 ~ "Stromal",
#                                                      `RNA_snn_res.0.2` == 10 ~ "Epithelial",
#                                                      `RNA_snn_res.0.2` == 11 ~ "Epithelial",
#                                                      `RNA_snn_res.0.2` == 12 ~ "Epithelial",
#                                                      `RNA_snn_res.0.2` == 13 ~ "Immune",
#                                                      `RNA_snn_res.0.2` == 14 ~ "Immune",
#                                                      `RNA_snn_res.0.2` == 15 ~ "Epithelial",
#                                                      `RNA_snn_res.0.2` == 16 ~ "Epithelial",
#                                                      `RNA_snn_res.0.2` == 17 ~ "Epithelial",
#                                                      `RNA_snn_res.0.2` == 18 ~ "Immune",
#                                                      `RNA_snn_res.0.2` == 19 ~ "Epithelial"))
combined_seurat_rpca$cluster_broad <- with(combined_seurat_rpca@meta.data,
                                           case_when(`RNA_snn_res.0.2` == 0 ~ "Immune",
                                                     `RNA_snn_res.0.2` == 1 ~ "Immune",
                                                     `RNA_snn_res.0.2` == 2 ~ "Epithelial",
                                                     `RNA_snn_res.0.2` == 3 ~ "Epithelial",
                                                     `RNA_snn_res.0.2` == 4 ~ "Immune",
                                                     `RNA_snn_res.0.2` == 5 ~ "Stromal",
                                                     `RNA_snn_res.0.2` == 6 ~ "Stromal",
                                                     `RNA_snn_res.0.2` == 7 ~ "Epithelial",
                                                     `RNA_snn_res.0.2` == 8 ~ "Epithelial",
                                                     `RNA_snn_res.0.2` == 9 ~ "Epithelial",
                                                     `RNA_snn_res.0.2` == 10 ~ "Stromal",
                                                     `RNA_snn_res.0.2` == 11 ~ "Immune",
                                                     `RNA_snn_res.0.2` == 12 ~ "Immune",
                                                     `RNA_snn_res.0.2` == 13 ~ "Epithelial",
                                                     `RNA_snn_res.0.2` == 14 ~ "Epithelial",
                                                     `RNA_snn_res.0.2` == 15 ~ "Epithelial",
                                                     `RNA_snn_res.0.2` == 16 ~ "Epithelial",
                                                     `RNA_snn_res.0.2` == 17 ~ "Immune",
                                                     `RNA_snn_res.0.2` == 18 ~ "Immune",
                                                     `RNA_snn_res.0.2` == 19 ~ "Immune",
                                                     `RNA_snn_res.0.2` == 20 ~ "Immune",
                                                     `RNA_snn_res.0.2` == 21 ~ "Stromal",
                                                     `RNA_snn_res.0.2` == 22 ~ "Immune",
                                                     `RNA_snn_res.0.2` == 23 ~ "Epithelial",
                                                     `RNA_snn_res.0.2` == 24 ~ "Epithelial"))
# combined_seurat_scanvi_ct$cluster_broad <- with(combined_seurat_scanvi_ct@meta.data,
#                                            case_when(`RNA_snn_res.0.2` == 0 ~ "Immune",
#                                                      `RNA_snn_res.0.2` == 1 ~ "Immune",
#                                                      `RNA_snn_res.0.2` == 2 ~ "Immune",
#                                                      `RNA_snn_res.0.2` == 3 ~ "Epithelial",
#                                                      `RNA_snn_res.0.2` == 4 ~ "Epithelial",
#                                                      `RNA_snn_res.0.2` == 5 ~ "Stromal",
#                                                      `RNA_snn_res.0.2` == 6 ~ "Stromal",
#                                                      `RNA_snn_res.0.2` == 7 ~ "Epithelial",
#                                                      `RNA_snn_res.0.2` == 8 ~ "Immune",
#                                                      `RNA_snn_res.0.2` == 9 ~ "Epithelial",
#                                                      `RNA_snn_res.0.2` == 10 ~ "Epithelial",
#                                                      `RNA_snn_res.0.2` == 11 ~ "Stromal",
#                                                      `RNA_snn_res.0.2` == 12 ~ "Immune",
#                                                      `RNA_snn_res.0.2` == 13 ~ "Epithelial",
#                                                      `RNA_snn_res.0.2` == 14 ~ "Epithelial",
#                                                      `RNA_snn_res.0.2` == 15 ~ "Epithelial",
#                                                      `RNA_snn_res.0.2` == 16 ~ "Epithelial",
#                                                      `RNA_snn_res.0.2` == 17 ~ "Epithelial",
#                                                      `RNA_snn_res.0.2` == 18 ~ "Epithelial"))
# 
# # 2. Confusion Matrix
# library(caret)
# harmony_rpca <- confusionMatrix(factor(combined_seurat_harmony$cluster_broad), 
#                                 factor(combined_seurat_rpca$cluster_broad))
# harmony_scanvi <- confusionMatrix(factor(combined_seurat_harmony$cluster_broad), 
#                                   factor(combined_seurat_scanvi_ct$cluster_broad))
# rpca_scanvi <- confusionMatrix(factor(combined_seurat_rpca$cluster_broad), 
#                                factor(combined_seurat_scanvi_ct$cluster_broad))
# 
# 
# # Plot the confusion matrix
# ggplot(as.data.frame(harmony_rpca$table), aes(x = Reference, y = Prediction)) +
#   geom_tile(aes(fill = Freq), color = "white") +
#   scale_fill_gradient(low = "white", high = "blue") +
#   geom_text(aes(label = Freq), vjust = 1) +
#   theme_minimal(base_size = 15) +
#   labs(title = "Confusion Matrix", x = "Harmony", y = "RPCA")
# ggplot(as.data.frame(harmony_scanvi$table), aes(x = Reference, y = Prediction)) +
#   geom_tile(aes(fill = Freq), color = "white") +
#   scale_fill_gradient(low = "white", high = "blue") +
#   geom_text(aes(label = Freq), vjust = 1) +
#   theme_minimal(base_size = 15) +
#   labs(title = "Confusion Matrix", x = "Harmony", y = "scANVI")
# ggplot(as.data.frame(rpca_scanvi$table), aes(x = Reference, y = Prediction)) +
#   geom_tile(aes(fill = Freq), color = "white") +
#   scale_fill_gradient(low = "white", high = "blue") +
#   geom_text(aes(label = Freq), vjust = 1) +
#   theme_minimal(base_size = 15) +
#   labs(title = "Confusion Matrix", x = "RPCA", y = "scANVI")
# 
# # 3. Proceeding with RPCA
# 
# # Comparing broad annotations
# p1 <- DimPlot_scCustom(combined_seurat_rpca,
#                        reduction = "umap.rpca",
#                        group.by = "RNA_snn_res.0.2",
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
# p5 <- FeaturePlot_scCustom(combined_seurat_rpca,
#                            reduction = "umap.rpca",
#                            features = "EPCAM") + labs(title = "EPCAM", x = "UMAP 1", y = "UMAP 2")
# p6 <- FeaturePlot_scCustom(combined_seurat_rpca,
#                        reduction = "umap.rpca",
#                        features = "PTPRC") + labs(title = "PTPRC", x = "UMAP 1", y = "UMAP 2")
# p7 <- FeaturePlot_scCustom(combined_seurat_rpca,
#                        reduction = "umap.rpca",
#                        features = "COL1A1") + labs(title = "COL1A1", x = "UMAP 1", y = "UMAP 2")
# p8 <- DimPlot_scCustom(combined_seurat_rpca,
#                        reduction = "umap.rpca",
#                        group.by="cluster_broad",
#                        raster = FALSE) + labs(title = "Final Compartments", x = "UMAP 1", y = "UMAP 2")
# combined <- cowplot::plot_grid(p1, p2, p3, p4, 
#                                p5, p6, p7, p8, ncol= 4, nrow = 2)
# 
# ggsave(plot = combined, filename = file.path(PATH, "results/umaps/rpca_compartments_all.png"),
#        width = 18, height = 8)
# imm_subset <- readRDS(file.path(PATH, "data/sc/imm_rpca_subset.rds"))
# epi_subset <- readRDS(file.path(PATH, "data/sc/epi_rpca_subset.rds"))
# strom_subset <- readRDS(file.path(PATH, "data/sc/strom_rpca_subset.rds"))
## Clustering each broad group
epi_subset <- subset(combined_seurat_rpca,
                     (cluster_broad == "Epithelial"))
epi_subset <- FindNeighbors(epi_subset, dims = 1:200, reduction = "integrated.rpca")
for (res in c(0.2, 0.4, 0.6, 0.8)) {
  epi_subset <- FindClusters(epi_subset, resolution = res)
}
# imm_subset <- subset(combined_seurat_rpca,
#                      (cluster_broad == "Immune"))
# imm_subset <- FindNeighbors(imm_subset, dims = 1:200, reduction = "integrated.rpca")
# for (res in c(0.2, 0.4, 0.6, 0.8)) {
#   imm_subset <- FindClusters(imm_subset, resolution = res)
# }
# strom_subset <- subset(combined_seurat_rpca,
#                        (cluster_broad == "Stromal"))
# strom_subset <- FindNeighbors(strom_subset, dims = 1:200, reduction = "integrated.rpca")
# for (res in c(0.2, 0.4, 0.6, 0.8)) {
#   strom_subset <- FindClusters(strom_subset, resolution = res)
# }

# imm_subset$celltypist_shortened <- celltypist_short(imm_subset$celltypist_pred)
# strom_subset$celltypist_shortened <- celltypist_short(strom_subset$celltypist_pred)
# epi_subset$celltypist_shortened <- celltypist_short(epi_subset$celltypist_pred)


#imm_subset <- RunUMAP(imm_subset, reduction = "integrated.rpca", dims = 1:200, reduction.name = "umap.rpca")
epi_subset <- RunUMAP(epi_subset, reduction = "integrated.rpca", dims = 1:200, reduction.name = "umap.rpca")
#strom_subset <- RunUMAP(strom_subset, reduction = "integrated.rpca", dims = 1:200, reduction.name = "umap.rpca")

## Celltypist shortened plots
# p1 <- DimPlot_scCustom(imm_subset,
#                        reduction = "umap.rpca",
#                        group.by = "celltypist_shortened",
#                        raster = FALSE) + labs(title = "Immune", x = "UMAP 1", y = "UMAP 2")
# ggsave(plot = p1, file.path(PATH, "results/umaps/imm_rpca_celltypist.png"), width = 10, height = 6)
# 
# p2 <- DimPlot_scCustom(strom_subset,
#                        reduction = "umap.rpca",
#                        group.by = "celltypist_shortened",
#                        raster = FALSE) + labs(title = "Stromal", x = "UMAP 1", y = "UMAP 2")
# ggsave(plot = p2, file.path(PATH, "results/umaps/strom_rpca_celltypist.png"), width = 10, height = 6)
# 
# p3 <- DimPlot_scCustom(epi_subset,
#                        reduction = "umap.rpca",
#                        group.by = "celltypist_shortened",
#                        raster = FALSE) + labs(title = "Epi", x = "UMAP 1", y = "UMAP 2")
# ggsave(plot = p3, file.path(PATH, "results/umaps/epi_rpca_celltypist.png"), width = 10, height = 6)

epi_subset <- JoinLayers(epi_subset)
saveRDS(epi_subset, file.path(PATH, "data/sc/epi_rpca_subset.rds"))
imm_subset <- readRDS(file.path(PATH, "data/sc/imm_rpca_subset.rds"))
imm_subset <- JoinLayers(imm_subset)
saveRDS(imm_subset, file.path(PATH, "data/sc/imm_rpca_subset.rds"))
#saveRDS(strom_subset, file.path(PATH, "data/sc/strom_rpca_subset.rds"))
print("Done with epi and  imm")
combined_seurat_rpca <- JoinLayers(combined_seurat_rpca)
saveRDS(combined_seurat_rpca, file.path(PATH, "data/sc/combined_seurat_rpca.rds"))

