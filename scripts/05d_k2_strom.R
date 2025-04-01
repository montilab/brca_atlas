library(Seurat)
library(tidyverse)
library(K2Taxonomer)
library(ggdendro)
library(Biobase)
library(hypeR)
library(AUCell)
#options(future.globals.maxSize = 100000 * 1024^2)
# print(unix::rlimit_all())
# unix::rlimit_as(100e12, 100e12)
# print(unix::rlimit_all())
source("cKmeansWrapperSqrt.R")

PATH <- file.path(Sys.getenv("MLAB"), "projects/brcameta/brca_atlas/")

strom_subset <- readRDS(file.path(PATH, "data/sc/strom_rpca_subset.rds"))
# Idents(strom_subset) <- strom_subset$RNA_snn_res.0.4
# 
# strom_subset$cluster_annot <- with(strom_subset@meta.data,
#                                       case_when(`RNA_snn_res.0.4` == 0 ~ "CAFs (COL11A1+)",
#                                                 `RNA_snn_res.0.4` == 1 ~ "CAFs (LAMP5)",
#                                                 `RNA_snn_res.0.4` == 2 ~ "Endo Imm.",
#                                                 `RNA_snn_res.0.4` == 3 ~ "Endo Ven.",
#                                                 `RNA_snn_res.0.4` == 4 ~ "VSMC",
#                                                 `RNA_snn_res.0.4` == 5 ~ "Endo Capil.",
#                                                 `RNA_snn_res.0.4` == 6 ~ "CAFs (PI16+)",
#                                                 `RNA_snn_res.0.4` == 7 ~ "Myofibroblast",
#                                                 `RNA_snn_res.0.4` == 8 ~ "CAFs (DPT+)",
#                                                 `RNA_snn_res.0.4` == 9 ~ "PCs (ECM)",
#                                                 `RNA_snn_res.0.4` == 10 ~ "Endo Arter.",
#                                                 `RNA_snn_res.0.4` == 11 ~ "CAFs (CA12+)",
#                                                 `RNA_snn_res.0.4` == 12 ~ "Endo Arter.",
#                                                 `RNA_snn_res.0.4` == 13 ~ "Endo Lymph.",
#                                                 `RNA_snn_res.0.4` == 14 ~ "Prolif Stromal",
#                                                 `RNA_snn_res.0.4` == 15 ~ "VSMC",
#                                                 `RNA_snn_res.0.4` == 16 ~ "CAFs (DPT+)"))
# saveRDS(strom_subset, file.path(PATH, "data/sc/strom_rpca_subset.rds"))

strom_subset$cluster_annot_num <- with(strom_subset@meta.data,
                                   case_when(`RNA_snn_res.0.4` == 12 ~ 10,
                                             `RNA_snn_res.0.4` == 15 ~ 4,
                                             `RNA_snn_res.0.4` == 16 ~ 8,
                                             .default = as.numeric(as.character(`RNA_snn_res.0.4`))))
strom_subset$cluster_annot_num <- as.factor(strom_subset$cluster_annot_num)
pan_cancer_strom_sigs <- readRDS(file.path(PATH, "data/signatures/stromal/ye_2024/pan_cancer_stromal_sigs.rds"))

RNGkind("L'Ecuyer-CMRG")
set.seed(1)
K2res <- K2preproc(object = t(strom_subset@reductions$integrated.rpca@cell.embeddings),
                   eMatDS = strom_subset@assays$RNA$data,
                   colData = strom_subset@meta.data,
                   variables = "batch",
                   cohorts="cluster_annot_num",
                   featMetric="F",
                   logCounts=TRUE,
                   clustFunc="cKmeansDownsampleSqrt",
                   useCors=20,
                   DGEmethod = "mast",
                   genesets = pan_cancer_strom_sigs,
                   ScoreGeneSetMethod = "AUCELL")

K2res <- K2tax(K2res)
K2res <- runDGEmods(K2res)
K2res <- runFISHERmods(K2res)
K2res <- runScoreGeneSets(K2res)
K2res <- runScoreGeneSets(K2res)
K2res <- runDSSEmods(K2res)
K2dashboard(K2res,  analysis_name = "K2Taxonomer_Strom")

# # Convert to expression set
# eSet <- ExpressionSet(assayData=as.matrix(t(strom_subset@reductions$integrated.rpca@cell.embeddings)))
# pData(eSet) <- as.data.frame(strom_subset@meta.data)
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
#                    covariates = "batch",
#                    logCounts=TRUE,
#                    nBoots=100,
#                    clustFunc=cKmeansWrapperSqrt,
#                    clustList=wrapperList)
# K2res <- K2tax(K2res)
# saveRDS(K2res, file.path(PATH, "data/k2_objects/k2_taxnomer_strom.rds"))
# # Get dendrogram from K2Taxonomer
# dendro <- K2dendro(K2res)
# 
# # Get dendrogram data
# dendro_plot <- ggdendrogram(dendro)
# ggsave(file.path(PATH, "results/k2/k2_dendro_strom.png"), plot = dendro_plot)
# 
# K2res <- readRDS(file.path(PATH, "data/k2_objects/k2_taxnomer_strom.rds"))
# # Annotation
# eSet <- ExpressionSet(assayData=as.matrix(strom_subset@assays$RNA$data))
# pData(eSet) <- as.data.frame(strom_subset@meta.data)
# 
# K2eSet(K2res) <- eSet
# K2meta(K2res)$covariates <- "batch"
# 
# K2res <- runDGEmods(K2res)
# ## Concatenate all results and generate table
# DGEtable <- getDGETable(K2res)
# head(DGEtable)
# saveRDS(K2res, file.path(PATH, "data/k2_objects/k2_taxnomer_strom.rds"))
# 
# # Enrichment
# K2res <- readRDS(file.path(PATH, "data/k2_objects/k2_taxnomer_strom.rds"))
# #strom_celltypes <- msigdb_download(species="Homo sapiens", category="H")
# strom_celltypes <- readRDS(file.path(PATH, "brca_atlas_validation/data/sigs/strom_markers.rds"))
# imm_celltypes <- readRDS(file.path(PATH, "brca_atlas_validation/data/sigs/imm_markers.rds"))
# wu_embo_strom <- readRDS(file.path(PATH, "data/signatures/stromal/wu_embo_2020/strom_sigs.rds"))
# strom_celltypes <- c(strom_celltypes, imm_celltypes[c("Pericyte", "VSMC")], wu_embo_strom)
# # Run gene set hyperenrichment
# K2res <- runGSEmods(K2res,
#                     genesets=strom_celltypes,
#                     qthresh=0.1)
# 
# K2_seurat <- CreateSeuratObject(counts = Biobase::exprs(K2eSet(K2res)),
#                                 meta.data = Biobase::pData(K2eSet(K2res)))
# 
# # AUCell
# K2_sce <- as.SingleCellExperiment(K2_seurat)
# cells_auc <- AUCell_run(K2_sce, strom_celltypes)
# strom_scores <- cells_auc@assays@data$AUC
# 
# strom_eSet <- ExpressionSet(assayData = as.matrix(strom_scores))
# pData(strom_eSet) <- pData(K2eSet(K2res))
# 
# K2gSet(K2res) <- strom_eSet
# 
# K2res <- runDSSEmods(K2res)
# K2dashboard(K2res, analysis_name = "K2Taxonomer_Strom_AUC")
# 
# # # VAM
# # K2_seurat@assays$RNA$data <- K2_seurat@assays$RNA$counts
# # K2_seurat <- split(K2_seurat, f = K2_seurat$batch)
# # K2_seurat <- FindVariableFeatures(K2_seurat)
# # K2_seurat <- JoinLayers(K2_seurat)
# # vam_geneset <- lapply(strom_celltypes, function(x) which(rownames(K2_seurat) %in% x))
# # vam_data <- VAM::vamForSeurat(K2_seurat, gene.set.collection = vam_geneset)
# # strom_scores <- vam_data@assays$VAMcdf$data
# # # Add Module Score
# # K2_seurat <- AddModuleScore(K2_seurat,
# #                             strom_celltypes,
# #                             slot="counts")
# # K2_seurat@meta.data <- K2_seurat@meta.data %>% data.table::setnames(old = paste0("Cluster",as.character(seq(1:length(strom_celltypes)))),
# #                                                                     new= names(strom_celltypes))
# # strom_scores <- K2_seurat@meta.data %>% dplyr::select(names(strom_celltypes))
# # strom_eSet <- ExpressionSet(assayData = t(as.matrix(strom_scores)))
# # pData(strom_eSet) <- pData(K2eSet(K2res))
# # 
# # K2gSet(K2res) <- strom_eSet
# # 
# # K2res <- runDSSEmods(K2res)
# # K2dashboard(K2res, analysis_name = "K2Taxonomer_Strom_AMS")