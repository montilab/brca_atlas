library(Seurat)
library(tidyverse)
library(K2Taxonomer)
library(ggdendro)
library(Biobase)
library(AUCell)
# source("cKmeansWrapperSqrt.R")
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

imm_subset <- readRDS(file.path(PATH, "data/sc/imm_rpca_sub_clean.rds"))

# imm_subset$cluster_annot_num <- with(imm_subset@meta.data,
#                                      case_when(`RNA_snn_res.0.6` == 28 ~ 7,
#                                                `RNA_snn_res.0.6` == 27 ~ 14,
#                                                `RNA_snn_res.0.6` == 24 ~ 4,
#                                                `RNA_snn_res.0.6` == 18 ~ 10,
#                                                .default = as.numeric(as.character(`RNA_snn_res.0.6`))))
# imm_subset$cluster_annot_num <- as.factor(imm_subset$cluster_annot_num)
immune_celltypes <- readRDS(file.path(PATH, "brca_atlas_validation/data/sigs/imm_markers.rds"))
boroni_sigs <- readRDS(file.path(PATH, "data/signatures/immune/boroni_2024/pan_cancer_myeloid_sigs_full.rds"))
vanderleun_sigs <- readRDS(file.path(PATH, "data/signatures/immune/vanderleun_sigs.rds"))
chu_sigs <- readRDS(file.path(PATH, "data/signatures/immune/chu_2023/chu_sigs.rds"))
all_imm_sigs <- c(immune_celltypes, boroni_sigs, vanderleun_sigs, chu_sigs)
all_imm_sigs <- all_imm_sigs[!duplicated(names(all_imm_sigs))]
saveRDS(all_imm_sigs, file.path(PATH, "data/signatures/immune/boroni_vanderleun_chu_hbca_sigs.rds"))

RNGkind("L'Ecuyer-CMRG")
set.seed(1)
K2res <- K2preproc(object = t(imm_subset@reductions$integrated.rpca@cell.embeddings),
                   eMatDS = imm_subset@assays$RNA$data,
                   colData = imm_subset@meta.data,
                   cohorts="cluster_annot_sub",
                   variables = "batch",
                   featMetric="F",
                   logCounts=TRUE,
                   clustFunc="cKmeansDownsampleSmallest",
                   nBoots=500,
                   useCors=28,
                   DGEmethod = "mast",
                   genesets = all_imm_sigs,
                   ScoreGeneSetMethod = "AUCELL")

K2res <- K2tax(K2res)
K2res <- runDGEmods(K2res)
K2res <- runFISHERmods(K2res)
K2res <- runScoreGeneSets(K2res)
K2res <- runDSSEmods(K2res)
K2dashboard(K2res,  analysis_name = "K2Taxonomer_Imm")
