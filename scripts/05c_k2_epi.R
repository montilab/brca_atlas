library(Seurat)
library(tidyverse)
library(K2Taxonomer)
library(ggdendro)
library(Biobase)
library(hypeR)
library(AUCell)
source("cKmeansWrapperSqrt.R")

PATH <- file.path(Sys.getenv("MLAB"), "projects/brcameta/brca_atlas/")

# epi_subset <- readRDS(file.path(PATH, "data/sc/epi_rpca_subset.rds"))
# epi_subset <- subset(epi_subset,
#                      subset = RNA_snn_res.0.2 %in% 0:14)
#   
# # Convert to expression set
# eSet <- ExpressionSet(assayData=as.matrix(t(epi_subset@reductions$integrated.rpca@cell.embeddings)))
# pData(eSet) <- as.data.frame(epi_subset@meta.data)
# 
# # Create clustList
# wrapperList <- list(
#   eMat=exprs(eSet),
#   labs=eSet$RNA_snn_res.0.2,
#   maxIter=10
# )
# 
# # Run K2Taxonomer
# K2res <- K2preproc(eSet,
#                    cohorts="RNA_snn_res.0.2",
#                    featMetric="F",
#                    logCounts=TRUE,
#                    nBoots=100,
#                    clustFunc=cKmeansWrapperSqrt,
#                    clustList=wrapperList)
# K2res <- K2tax(K2res)
# saveRDS(K2res, file.path(PATH, "data/k2_objects/k2_taxnomer_epi.rds"))
# # Get dendrogram from K2Taxonomer
# dendro <- K2dendro(K2res)
# 
# # Get dendrogram data
# dendro_plot <- ggdendrogram(dendro)
# ggsave(file.path(PATH, "results/k2/k2_dendro_epi.png"), plot = dendro_plot)

# Annotation
K2res <- readRDS(file.path(PATH, "data/k2_objects/k2_taxnomer_epi.rds"))
eSet <- ExpressionSet(assayData=as.matrix(epi_subset@assays$RNA$data))
pData(eSet) <- as.data.frame(epi_subset@meta.data)
K2eSet(K2res) <- eSet

K2meta(K2res)$covariates <- "batch"

K2res <- runDGEmods(K2res)
## Concatenate all results and generate table
DGEtable <- getDGETable(K2res)
head(DGEtable)
saveRDS(K2res, file.path(PATH, "data/k2_objects/k2_taxnomer_epi.rds"))

# Enrichment
K2res <- readRDS(file.path(PATH, "data/k2_objects/k2_taxnomer_epi.rds"))
# epi_celltypes <- msigdb_download(species="Homo sapiens", category="H")
# epi_celltypes <- readRDS(file.path(PATH, "brca_atlas_validation/data/sigs/epi_markers.rds"))
# epi_barkley <- readRDS(file.path(PATH, "data/signatures/cancer_epithelial/barkley_natgen_2023/barkley_sigs.rds"))
# epi_gavish <- readRDS(file.path(PATH, "data/signatures/cancer_epithelial/gavish_nature_2023/gavish_sigs.rds"))
# epi_xu <- readRDS(file.path(PATH, "data/signatures/cancer_epithelial/xu_cellreports_2024/xu_sigs.rds"))
# epi_all <- c(epi_celltypes, epi_barkley, epi_gavish, epi_xu)
epi_all <- readRDS(file.path(PATH, "brca_atlas_validation/data/sigs/all_markers.rds"))
epi_all <- epi_all[1:11]
# Run gene set hyperenrichment
K2res <- runGSEmods(K2res,
                    genesets=epi_all,
                    qthresh=0.1)
#counts here is actually normalized data
K2_seurat <- CreateSeuratObject(counts = Biobase::exprs(K2eSet(K2res)),
                                meta.data = Biobase::pData(K2eSet(K2res)))
# AUC
K2_sce <- as.SingleCellExperiment(K2_seurat)
cells_auc <- AUCell_run(K2_sce, epi_all)
epi_scores <- cells_auc@assays@data$AUC
# # VAM Scores
# K2_seurat@assays$RNA$data <- K2_seurat@assays$RNA$counts
# K2_seurat <- split(K2_seurat, f = K2_seurat$batch)
# K2_seurat <- FindVariableFeatures(K2_seurat)
# K2_seurat <- JoinLayers(K2_seurat)
# vam_geneset <- lapply(epi_all, function(x) which(rownames(K2_seurat) %in% x))
# vam_data <- VAM::vamForSeurat(K2_seurat, gene.set.collection = vam_geneset)
# epi_scores <- vam_data@assays$VAMcdf$data
# # Add Module Score
# K2_seurat <- AddModuleScore(K2_seurat,
#                             epi_all,
#                             slot = "counts")
# K2_seurat@meta.data <- K2_seurat@meta.data %>% data.table::setnames(old = paste0("Cluster",as.character(seq(1:length(epi_all)))),
#                                                                     new= names(epi_all))
#epi_scores <- K2_seurat@meta.data %>% dplyr::select(names(epi_all))
epi_eSet <- ExpressionSet(assayData = as.matrix(epi_scores))
pData(epi_eSet) <- pData(K2eSet(K2res))

K2gSet(K2res) <- epi_eSet

K2res <- runDSSEmods(K2res)
K2dashboard(K2res, analysis_name = "K2Taxonomer_Epi_subset_AUC")

# # Add Module Score
# K2_seurat <- AddModuleScore(K2_seurat,
#                             epi_all,
#                             slot = "counts")
# K2_seurat@meta.data <- K2_seurat@meta.data %>% data.table::setnames(old = paste0("Cluster",as.character(seq(1:length(epi_all)))),
#                                                                     new= names(epi_all))
# epi_scores <- K2_seurat@meta.data %>% dplyr::select(names(epi_all))
# epi_eSet <- ExpressionSet(assayData = t(as.matrix(epi_scores)))
# pData(epi_eSet) <- pData(K2eSet(K2res))
# 
# K2gSet(K2res) <- epi_eSet
# 
# K2res <- runDSSEmods(K2res)
# K2dashboard(K2res, analysis_name = "K2Taxonomer_Epi_AMS")
