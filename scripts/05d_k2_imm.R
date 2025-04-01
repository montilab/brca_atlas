library(Seurat)
library(tidyverse)
library(K2Taxonomer)
library(ggdendro)
library(Biobase)
library(AUCell)
source("cKmeansWrapperSqrt.R")
options(future.globals.maxSize = 8000 * 1024^2)
PATH <- file.path(Sys.getenv("MLAB"), "projects/brcameta/brca_atlas/")

# imm_rpca_subset <- readRDS(file.path(PATH, "data/sc/imm_rpca_subset.rds"))
# imm_rpca_subset$cluster_annot <- with(imm_rpca_subset@meta.data,
#                                       case_when(`RNA_snn_res.0.6` == 0 ~ "Mac Lipo",
#                                                 `RNA_snn_res.0.6` == 1 ~ "CD8 Tem",
#                                                 `RNA_snn_res.0.6` == 2 ~ "CD4 Tem",
#                                                 `RNA_snn_res.0.6` == 3 ~ "NK",
#                                                 `RNA_snn_res.0.6` == 4 ~ "B Cell",
#                                                 `RNA_snn_res.0.6` == 5 ~ "CD4 Treg",
#                                                 `RNA_snn_res.0.6` == 6 ~ "CD4 Naive",
#                                                 `RNA_snn_res.0.6` == 7 ~ "Plasma IgG",
#                                                 `RNA_snn_res.0.6` == 8 ~ "Mac Col.",
#                                                 `RNA_snn_res.0.6` == 9 ~ "CD8 Isg",
#                                                 `RNA_snn_res.0.6` == 10 ~ "CD8 Tex",
#                                                 `RNA_snn_res.0.6` == 11 ~ "cDC2",
#                                                 `RNA_snn_res.0.6` == 12 ~ "CD4 Tfh",
#                                                 `RNA_snn_res.0.6` == 13 ~ "Other",
#                                                 `RNA_snn_res.0.6` == 14 ~ "T Prolif.",
#                                                 `RNA_snn_res.0.6` == 15 ~ "Other",
#                                                 `RNA_snn_res.0.6` == 16 ~ "Mac Prolif.",
#                                                 `RNA_snn_res.0.6` == 17 ~ "Mast",
#                                                 `RNA_snn_res.0.6` == 18 ~ "CD8 Tex",
#                                                 `RNA_snn_res.0.6` == 19 ~ "pDC",
#                                                 `RNA_snn_res.0.6` == 20 ~ "Other",
#                                                 `RNA_snn_res.0.6` == 21 ~ "mDC",
#                                                 `RNA_snn_res.0.6` == 22 ~ "Other",
#                                                 `RNA_snn_res.0.6` == 23 ~ "cDC1",
#                                                 `RNA_snn_res.0.6` == 24 ~ "B Cell",
#                                                 `RNA_snn_res.0.6` == 25 ~ "CD4 Th",
#                                                 `RNA_snn_res.0.6` == 26 ~ "Other",
#                                                 `RNA_snn_res.0.6` == 27 ~ "T Prolif.",
#                                                 `RNA_snn_res.0.6` == 28 ~ "Plasma IgG"))
# 
# imm_subset <- subset(imm_rpca_subset,
#                          (cluster_annot != "Other"))
# imm_subset$`RNA_snn_res.0.6` <- droplevels(imm_subset$`RNA_snn_res.0.6`)
# saveRDS(imm_subset, file.path(PATH, "data/sc/imm_rpca_clean.rds"))

imm_subset <- readRDS(file.path(PATH, "data/sc/imm_rpca_clean.rds"))

imm_subset$cluster_annot_num <- with(imm_subset@meta.data,
                                     case_when(`RNA_snn_res.0.6` == 28 ~ 7,
                                               `RNA_snn_res.0.6` == 27 ~ 14,
                                               `RNA_snn_res.0.6` == 24 ~ 4,
                                               `RNA_snn_res.0.6` == 18 ~ 10,
                                               .default = as.numeric(as.character(`RNA_snn_res.0.6`))))
imm_subset$cluster_annot_num <- as.factor(imm_subset$cluster_annot_num)
all_imm_sigs <- readRDS(file.path(PATH, "data/signatures/immune/boroni_vanderleun_chu_hbca_sigs.rds"))

RNGkind("L'Ecuyer-CMRG")
set.seed(1)
K2res <- K2preproc(object = t(imm_subset@reductions$integrated.rpca@cell.embeddings),
                   eMatDS = imm_subset@assays$RNA$data,
                   colData = imm_subset@meta.data,
                   cohorts="cluster_annot_num",
                   variables = "batch",
                   featMetric="F",
                   logCounts=TRUE,
                   clustFunc="cKmeansDownsampleSqrt",
                   useCors=20,
                   DGEmethod = "mast",
                   genesets = all_imm_sigs,
                   ScoreGeneSetMethod = "AUCELL")

K2res <- K2tax(K2res)
K2res <- runDGEmods(K2res)
K2res <- runFISHERmods(K2res)
K2res <- runScoreGeneSets(K2res)
K2res <- runScoreGeneSets(K2res)
K2res <- runDSSEmods(K2res)
K2dashboard(K2res,  analysis_name = "K2Taxonomer_Imm")

# # Convert to expression set
# eSet <- ExpressionSet(assayData=as.matrix(t(imm_subset@reductions$integrated.rpca@cell.embeddings)))
# pData(eSet) <- as.data.frame(imm_subset@meta.data)
# 
# # Create clustList
# wrapperList <- list(
#   eMat=exprs(eSet),
#   labs=eSet$cluster_annot_num,
#   maxIter=10
# )
# 
# # Run K2Taxonomer
# K2res <- K2preproc(eSet,
#                    cohorts="cluster_annot_num",
#                    featMetric="F",
#                    logCounts=TRUE,
#                    nBoots=100,
#                    clustFunc=cKmeansWrapperSqrt,
#                    clustList=wrapperList)
# K2res <- K2tax(K2res)
# saveRDS(K2res, file.path(PATH, "data/k2_objects/k2_taxnomer_imm.rds"))
# ## Get dendrogram from K2Taxonomer
# dendro <- K2dendro(K2res)
# 
# ## Get dendrogram data
# dendro_plot <- ggdendrogram(dendro)
# ggsave(file.path(PATH, "results/k2/k2_dendro_imm.png"), plot = dendro_plot)
# print("Done with dendrogram")
# 
# # Annotation
# K2res <- readRDS(file.path(PATH, "data/k2_objects/k2_taxnomer_imm.rds"))
# eSet <- ExpressionSet(assayData=as.matrix(imm_subset@assays$RNA$data))
# pData(eSet) <- as.data.frame(imm_subset@meta.data)
# K2eSet(K2res) <- eSet
# 
# K2meta(K2res)$covariates <- "batch"
# 
# K2res <- runDGEmods(K2res)
# ## Concatenate all results and generate table
# DGEtable <- getDGETable(K2res)
# head(DGEtable)
# saveRDS(K2res, file.path(PATH, "data/k2_objects/k2_taxnomer_imm.rds"))
# print("Done with annotation")
# 
# # Enrichment
# K2res <- readRDS(file.path(PATH, "data/k2_objects/k2_taxnomer_imm.rds"))
# immune_celltypes <- readRDS(file.path(PATH, "brca_atlas_validation/data/sigs/imm_markers.rds"))
# boroni_sigs <- readRDS(file.path(PATH, "data/signatures/immune/boroni_2024/pan_cancer_myeloid_sigs.rds"))
# vanderleun_sigs <- readRDS(file.path(PATH, "data/signatures/immune/vanderleun_sigs.rds"))
# chu_sigs <- readRDS(file.path(PATH, "data/signatures/immune/chu_2023/chu_sigs.rds"))
# all_imm_sigs <- c(immune_celltypes, boroni_sigs, vanderleun_sigs, chu_sigs)
# all_imm_sigs <- all_imm_sigs[!duplicated(names(all_imm_sigs))]
# saveRDS(all_imm_sigs, file.path(PATH, "data/signatures/immune/boroni_vanderleun_chu_hbca_sigs.rds"))
# # all_imm_sigs <- all_imm_sigs[as.numeric(lapply(all_imm_sigs, length)) > 0]
# # all_imm_sigs <- lapply(all_imm_sigs, function(x) x[x %in% rownames(imm_subset)])
# 
# ## Run gene set hyperenrichment
# K2res <- runGSEmods(K2res,
#                     genesets=all_imm_sigs,
#                     qthresh=0.1)
# K2_seurat <- CreateSeuratObject(counts = Biobase::exprs(K2eSet(K2res)),
#                                meta.data = Biobase::pData(K2eSet(K2res)))
# 
# # AUC
# K2_sce <- as.SingleCellExperiment(K2_seurat)
# cells_auc <- AUCell_run(K2_sce, all_imm_sigs)
# immune_scores <- cells_auc@assays@data$AUC
# # # VAM
# # K2_seurat@assays$RNA$data <- K2_seurat@assays$RNA$counts
# # K2_seurat <- split(K2_seurat, f = K2_seurat$batch)
# # K2_seurat <- FindVariableFeatures(K2_seurat)
# # K2_seurat <- JoinLayers(K2_seurat)
# # vam_geneset <- lapply(immune_celltypes, function(x) which(rownames(K2_seurat) %in% x))
# # vam_data <- VAM::vamForSeurat(K2_seurat, gene.set.collection = vam_geneset)
# # immune_scores <- vam_data@assays$VAMcdf$data
# 
# immune_eset <- ExpressionSet(assayData = as.matrix(immune_scores))
# pData(immune_eset) <- pData(K2eSet(K2res))
# 
# K2gSet(K2res) <- immune_eset
# 
# K2res <- runDSSEmods(K2res)
# K2dashboard(K2res, analysis_name = "K2Taxonomer_Imm_AUC")
# 
# # # Add Module Score
# # K2_seurat <- AddModuleScore(K2_seurat,
# #                             immune_celltypes,
# #                             slot = "counts")
# # K2_seurat@meta.data <- K2_seurat@meta.data %>% data.table::setnames(old = paste0("Cluster",as.character(seq(1:length(immune_celltypes)))),
# #                                                                     new= names(immune_celltypes))
# # immune_scores <- K2_seurat@meta.data %>% dplyr::select(names(immune_celltypes))
# # immune_eset <- ExpressionSet(assayData = t(as.matrix(immune_scores)))
# # pData(immune_eset) <- pData(K2eSet(K2res))
# # 
# # K2gSet(K2res) <- immune_eset
# # 
# # K2res <- runDSSEmods(K2res)
# # K2dashboard(K2res, analysis_name = "K2Taxonomer_Imm_AMS")