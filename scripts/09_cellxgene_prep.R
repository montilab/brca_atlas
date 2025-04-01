library(Seurat)
library(tidyverse)
library(sceasy)
# library(biomaRt)
# library(AnnotationDbi)
# library(org.Hs.eg.db)
library(reticulate)
reticulate::use_condaenv("r-sceasy")
options(Seurat.object.assay.version = "v5")
# ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
PATH <- file.path(Sys.getenv("MLAB"), "projects/brcameta/brca_atlas/")
DATA_PATH <- file.path(Sys.getenv("CBM"), "otherStudies/scRNAseq/brca")

# Using Pal & Wang mapping tables
pal_feats <- read.delim(file.path(DATA_PATH, "pal_2021/GSE161529_features.tsv"), header = FALSE)[,1:2]
colnames(pal_feats) <- c("ensembl", "hgnc")
wang_feats <- read.delim(file.path(DATA_PATH, "wang_2024/TNBC/GSM8064219_BC17086-12_Tumor_features.tsv.gz"), header = FALSE)[,1:2]
colnames(wang_feats) <- c("ensembl", "hgnc")
wang_feats$ensembl <- str_remove_all(wang_feats$ensembl, pattern = "\\.\\d+$")

merged_ensembl_hgnc <- dplyr::full_join(pal_feats, wang_feats, by = "ensembl")

cellxgene_metadata_columns <- c("organism_ontology_term_id",
                                "tissue_ontology_term_id",
                                "tissue_type",
                                "assay_ontology_term_id",
                                "disease_ontology_term_id",
                                "cell_type_ontology_term_id",
                                "self_reported_ethnicity_ontology_term_id",
                                "development_stage_ontology_term_id",
                                "sex_ontology_term_id",
                                "donor_id",
                                "suspension_type",
                                "age",
                                "grade",
                                "author_cell_type",
                                "batch")

# # 1. Imm Subset
# imm_rpca <- readRDS(file.path(PATH, "data/sc/imm_rpca_clean.rds"))
#
# # Filtering gene symbols to matchable ensemblids
# hgnc_ensembl <- data.frame(hgnc = rownames(imm_rpca))
# hgnc_ensembl$ensembl <- with(hgnc_ensembl, case_when(hgnc %in% merged_ensembl_hgnc$hgnc.x ~ merged_ensembl_hgnc$ensembl[match(hgnc, merged_ensembl_hgnc$hgnc.x)],
#                                                      hgnc %in% merged_ensembl_hgnc$hgnc.y ~ merged_ensembl_hgnc$ensembl[match(hgnc, merged_ensembl_hgnc$hgnc.y)],
#                                                      .default = NA))
# hgnc_ensembl_filtered <- hgnc_ensembl[!is.na(hgnc_ensembl$ensembl),]
# hgnc_ensembl_filtered <- hgnc_ensembl_filtered[!duplicated(hgnc_ensembl_filtered$ensembl),]
# rownames_filter <- rownames(imm_rpca) %in% hgnc_ensembl_filtered$hgnc
#
# # Subset and replace rownames
# imm_rpca <- imm_rpca[rownames_filter,]
# counts_layer <- imm_rpca@assays$RNA$counts
# rownames(counts_layer) <- hgnc_ensembl_filtered$ensembl[match(rownames(counts_layer), hgnc_ensembl_filtered$hgnc)]
# data_layer <- imm_rpca@assays$RNA$data
# rownames(data_layer) <- hgnc_ensembl_filtered$ensembl[match(rownames(data_layer), hgnc_ensembl_filtered$hgnc)]
# scale_data_layer <- imm_rpca@assays$RNA$scale.data
# rownames(scale_data_layer) <- hgnc_ensembl_filtered$ensembl[match(rownames(scale_data_layer), hgnc_ensembl_filtered$hgnc)]
# imm_rpca_new <- CreateSeuratObject(counts = counts_layer,
#                                      meta.data = imm_rpca@meta.data,
#                                      project = "cellxgene")
# imm_rpca_new@assays$RNA$data <- data_layer
# imm_rpca_new@assays$RNA$scale.data <- scale_data_layer
# for (reduction in names(imm_rpca@reductions)) {
#   imm_rpca_new[[reduction]] <- imm_rpca[[reduction]]
# }
#
# # Add metadata information
# imm_rpca_new$organism_ontology_term_id <- "NCBITaxon:9606"
# imm_rpca_new$tissue_ontology_term_id <- "UBERON:0000310"
# imm_rpca_new$tissue_type <- "tissue"
# imm_rpca_new$assay_ontology_term_id <- with(imm_rpca_new@meta.data,
#                                              case_when(batch == "wu_natgen_2021" ~ "EFO:0009899",
#                                                        batch == "pal_2021" ~ "EFO:0030003",
#                                                        batch == "qian_2020" ~ "EFO:0009900",
#                                                        batch == "gao_2021" ~ "EFO:0009899",
#                                                        batch == "tietscher_2023" ~ "EFO:0009922",
#                                                        batch == "liu_2023" ~ "EFO:0009922",
#                                                        batch == "wang_2024" ~ "EFO:0009899",
#                                                        batch == "bassez_2021" ~ "EFO:0030004",
#                                                        donor %in% c("wu_natgen_CID4040", "wu_natgen_CID3838", "wu_natgen_CID4040", "BC393_Tumor", "BC392_Tumor") ~ "EFO:0009922"))
# imm_rpca_new$disease_ontology_term_id <- "MONDO:0007254"
# # Match Cell type Ontologies to annotations found in dataset
# imm_rpca_new$cell_type_ontology_term_id <- with(imm_rpca_new@meta.data,
#                                       case_when(`RNA_snn_res.0.6` == 0 ~ "CL:0000235",
#                                                 `RNA_snn_res.0.6` == 1 ~ "CL:0000913",
#                                                 `RNA_snn_res.0.6` == 2 ~ "CL:0000905",
#                                                 `RNA_snn_res.0.6` == 3 ~ "CL:0000623",
#                                                 `RNA_snn_res.0.6` == 4 ~ "CL:0000787",
#                                                 `RNA_snn_res.0.6` == 5 ~ "CL:0001047",
#                                                 `RNA_snn_res.0.6` == 6 ~ "CL:0000898",
#                                                 `RNA_snn_res.0.6` == 7 ~ "CL:0000985",
#                                                 `RNA_snn_res.0.6` == 8 ~ "CL:0000235",
#                                                 `RNA_snn_res.0.6` == 9 ~ "CL:0000795",
#                                                 `RNA_snn_res.0.6` == 10 ~ "CL:0011025",
#                                                 `RNA_snn_res.0.6` == 11 ~ "CL:0000990",
#                                                 `RNA_snn_res.0.6` == 12 ~ "CL:0002038",
#                                                 `RNA_snn_res.0.6` == 14 ~ "CL:4033069",
#                                                 `RNA_snn_res.0.6` == 16 ~ "CL:4033076",
#                                                 `RNA_snn_res.0.6` == 17 ~ "CL:0000097",
#                                                 `RNA_snn_res.0.6` == 18 ~ "CL:0011025",
#                                                 `RNA_snn_res.0.6` == 19 ~ "CL:0000784",
#                                                 `RNA_snn_res.0.6` == 21 ~ "CL:0000782",
#                                                 `RNA_snn_res.0.6` == 23 ~ "CL:0000990",
#                                                 `RNA_snn_res.0.6` == 24 ~ "CL:0000787",
#                                                 `RNA_snn_res.0.6` == 25 ~ "CL:0000492",
#                                                 `RNA_snn_res.0.6` == 27 ~ "CL:4033069",
#                                                 `RNA_snn_res.0.6` == 28 ~ "CL:0000985"))
# imm_rpca_new$self_reported_ethnicity_ontology_term_id <- "unknown"
# # Cell x gene will help you translate numerical age values to ontology ids
# imm_rpca_new$development_stage_ontology_term_id <- imm_rpca_new$age
# imm_rpca_new$sex_ontology_term_id <- "PATO:0000383"
# imm_rpca_new$donor_id <- imm_rpca_new$donor
# imm_rpca_new$suspension_type <- "cell"
# imm_rpca_new$author_cell_type <- imm_rpca_new$cluster_annot
# imm_rpca_new@reductions$pca <- NULL
# imm_rpca_new@reductions$umap.pca <- NULL
# imm_rpca_new@reductions$integrated.rpca <- NULL
# imm_rpca_new[["RNA"]] <- as(imm_rpca_new[["RNA"]], "Assay")
# imm_rpca_new@meta.data <- imm_rpca_new@meta.data %>% dplyr::select(all_of(cellxgene_metadata_columns))
# saveRDS(imm_rpca_new, file.path(PATH, "data/sc/cellxgene/imm.rds"))
# imm_rpca_new <- readRDS(file.path(PATH, "data/sc/cellxgene/imm.rds"))
# sceasy::convertFormat(imm_rpca_new, from="seurat", to="anndata",
#                       main_layer = "counts", transfer_layers="data",
#                       drop_single_values = FALSE,
#                       outFile=file.path(PATH, "data/sc/cellxgene/imm.h5ad"))
#
# # 2. Strom Subset
# strom_rpca <- readRDS(file.path(PATH, "data/sc/strom_rpca_subset.rds"))
#
# # Filtering gene symbols to matchable ensemblids
# hgnc_ensembl <- data.frame(hgnc = rownames(strom_rpca))
# hgnc_ensembl$ensembl <- with(hgnc_ensembl, case_when(hgnc %in% merged_ensembl_hgnc$hgnc.x ~ merged_ensembl_hgnc$ensembl[match(hgnc, merged_ensembl_hgnc$hgnc.x)],
#                                                      hgnc %in% merged_ensembl_hgnc$hgnc.y ~ merged_ensembl_hgnc$ensembl[match(hgnc, merged_ensembl_hgnc$hgnc.y)],
#                                                      .default = NA))
# hgnc_ensembl_filtered <- hgnc_ensembl[!is.na(hgnc_ensembl$ensembl),]
# hgnc_ensembl_filtered <- hgnc_ensembl_filtered[!duplicated(hgnc_ensembl_filtered$ensembl),]
# rownames_filter <- rownames(strom_rpca) %in% hgnc_ensembl_filtered$hgnc
#
# # Subset and replace rownames
# strom_rpca <- strom_rpca[rownames_filter,]
# counts_layer <- strom_rpca@assays$RNA$counts
# rownames(counts_layer) <- hgnc_ensembl_filtered$ensembl[match(rownames(counts_layer), hgnc_ensembl_filtered$hgnc)]
# data_layer <- strom_rpca@assays$RNA$data
# rownames(data_layer) <- hgnc_ensembl_filtered$ensembl[match(rownames(data_layer), hgnc_ensembl_filtered$hgnc)]
# scale_data_layer <- strom_rpca@assays$RNA$scale.data
# rownames(scale_data_layer) <- hgnc_ensembl_filtered$ensembl[match(rownames(scale_data_layer), hgnc_ensembl_filtered$hgnc)]
# strom_rpca_new <- CreateSeuratObject(counts = counts_layer,
#                                      meta.data = strom_rpca@meta.data,
#                                      project = "cellxgene")
# strom_rpca_new@assays$RNA$data <- data_layer
# strom_rpca_new@assays$RNA$scale.data <- scale_data_layer
# for (reduction in names(strom_rpca@reductions)) {
#   strom_rpca_new[[reduction]] <- strom_rpca[[reduction]]
# }
# # Add metadata information
# strom_rpca_new$organism_ontology_term_id <- "NCBITaxon:9606"
# strom_rpca_new$tissue_ontology_term_id <- "UBERON:0000310"
# strom_rpca_new$tissue_type <- "tissue"
# strom_rpca_new$assay_ontology_term_id <- with(strom_rpca_new@meta.data,
#                                         case_when(batch == "wu_natgen_2021" ~ "EFO:0009899",
#                                                   batch == "pal_2021" ~ "EFO:0030003",
#                                                   batch == "qian_2020" ~ "EFO:0009900",
#                                                   batch == "gao_2021" ~ "EFO:0009899",
#                                                   batch == "tietscher_2023" ~ "EFO:0009922",
#                                                   batch == "liu_2023" ~ "EFO:0009922",
#                                                   batch == "wang_2024" ~ "EFO:0009899",
#                                                   batch == "bassez_2021" ~ "EFO:0030004",
#                                                   donor %in% c("wu_natgen_CID4040", "wu_natgen_CID3838", "wu_natgen_CID4040", "BC393_Tumor", "BC392_Tumor") ~ "EFO:0009922"))
# strom_rpca_new$disease_ontology_term_id <- "MONDO:0007254"
# # Match Cell type Ontologies to annotations found in dataset
# strom_rpca_new$cell_type_ontology_term_id <- with(strom_rpca_new@meta.data,
#                                             case_when(`RNA_snn_res.0.4` == 0 ~ "CL:0000057",
#                                                       `RNA_snn_res.0.4` == 1 ~ "CL:0000057",
#                                                       `RNA_snn_res.0.4` == 2 ~ "CL:0000115",
#                                                       `RNA_snn_res.0.4` == 3 ~ "CL:0002543",
#                                                       `RNA_snn_res.0.4` == 4 ~ "CL:0000359",
#                                                       `RNA_snn_res.0.4` == 5 ~ "CL:0002139",
#                                                       `RNA_snn_res.0.4` == 6 ~ "CL:0000057",
#                                                       `RNA_snn_res.0.4` == 7 ~ "CL:0000186",
#                                                       `RNA_snn_res.0.4` == 8 ~ "CL:0000057",
#                                                       `RNA_snn_res.0.4` == 9 ~ "CL:0000669",
#                                                       `RNA_snn_res.0.4` == 10 ~ "CL:1000413",
#                                                       `RNA_snn_res.0.4` == 11 ~ "CL:0000057",
#                                                       `RNA_snn_res.0.4` == 12 ~ "CL:0002144",
#                                                       `RNA_snn_res.0.4` == 13 ~ "CL:0000115",
#                                                       `RNA_snn_res.0.4` == 14 ~ "CL:4047001",
#                                                       `RNA_snn_res.0.4` == 15 ~ "CL:0000359",
#                                                       `RNA_snn_res.0.4` == 16 ~ "CL:0000057"))
# strom_rpca_new$self_reported_ethnicity_ontology_term_id <- "unknown"
# # Cell x gene will help you translate numerical age values to ontology ids
# strom_rpca_new$development_stage_ontology_term_id <- strom_rpca_new$age
# strom_rpca_new$sex_ontology_term_id <- "PATO:0000383"
# strom_rpca_new$donor_id <- strom_rpca_new$donor
# strom_rpca_new$suspension_type <- "cell"
# strom_rpca_new$author_cell_type <- strom_rpca_new$cluster_annot
# strom_rpca_new@reductions$pca <- NULL
# strom_rpca_new@reductions$umap.pca <- NULL
# strom_rpca_new@reductions$integrated.rpca <- NULL
# strom_rpca_new[["RNA"]] <- as(strom_rpca_new[["RNA"]], "Assay")
# strom_rpca_new@meta.data <- strom_rpca_new@meta.data %>% dplyr::select(all_of(cellxgene_metadata_columns))
# saveRDS(strom_rpca_new, file.path(PATH, "data/sc/cellxgene/strom.rds"))
# sceasy::convertFormat(strom_rpca_new, from="seurat", to="anndata",
#                       main_layer = "counts", transfer_layers="data",
#                       drop_single_values = FALSE,
#                       outFile=file.path(PATH, "data/sc/cellxgene/strom.h5ad"))
#
# 3. Epi Subset
epi_rpca <- readRDS(file.path(PATH, "data/sc/epi_rpca_subset.rds"))

# Filtering gene symbols to matchable ensemblids
hgnc_ensembl <- data.frame(hgnc = rownames(epi_rpca))
hgnc_ensembl$ensembl <- with(hgnc_ensembl, case_when(hgnc %in% merged_ensembl_hgnc$hgnc.x ~ merged_ensembl_hgnc$ensembl[match(hgnc, merged_ensembl_hgnc$hgnc.x)],
                                                     hgnc %in% merged_ensembl_hgnc$hgnc.y ~ merged_ensembl_hgnc$ensembl[match(hgnc, merged_ensembl_hgnc$hgnc.y)],
                                                     .default = NA))
hgnc_ensembl_filtered <- hgnc_ensembl[!is.na(hgnc_ensembl$ensembl),]
hgnc_ensembl_filtered <- hgnc_ensembl_filtered[!duplicated(hgnc_ensembl_filtered$ensembl),]
rownames_filter <- rownames(epi_rpca) %in% hgnc_ensembl_filtered$hgnc

# Subset and replace rownames
epi_rpca <- epi_rpca[rownames_filter,]
counts_layer <- epi_rpca@assays$RNA$counts
rownames(counts_layer) <- hgnc_ensembl_filtered$ensembl[match(rownames(counts_layer), hgnc_ensembl_filtered$hgnc)]
data_layer <- epi_rpca@assays$RNA$data
rownames(data_layer) <- hgnc_ensembl_filtered$ensembl[match(rownames(data_layer), hgnc_ensembl_filtered$hgnc)]
scale_data_layer <- epi_rpca@assays$RNA$scale.data
rownames(scale_data_layer) <- hgnc_ensembl_filtered$ensembl[match(rownames(scale_data_layer), hgnc_ensembl_filtered$hgnc)]
epi_rpca_new <- CreateSeuratObject(counts = counts_layer,
                                     meta.data = epi_rpca@meta.data,
                                     project = "cellxgene")
epi_rpca_new@assays$RNA$data <- data_layer
epi_rpca_new@assays$RNA$scale.data <- scale_data_layer
for (reduction in names(epi_rpca@reductions)) {
  epi_rpca_new[[reduction]] <- epi_rpca[[reduction]]
}

# Add metadata information
epi_rpca_new$organism_ontology_term_id <- "NCBITaxon:9606"
epi_rpca_new$tissue_ontology_term_id <- "UBERON:0000310"
epi_rpca_new$tissue_type <- "tissue"
epi_rpca_new$assay_ontology_term_id <- with(epi_rpca_new@meta.data,
                                        case_when(batch == "wu_natgen_2021" ~ "EFO:0009899",
                                                  batch == "pal_2021" ~ "EFO:0030003",
                                                  batch == "qian_2020" ~ "EFO:0009900",
                                                  batch == "gao_2021" ~ "EFO:0009899",
                                                  batch == "tietscher_2023" ~ "EFO:0009922",
                                                  batch == "liu_2023" ~ "EFO:0009922",
                                                  batch == "wang_2024" ~ "EFO:0009899",
                                                  batch == "bassez_2021" ~ "EFO:0030004",
                                                  donor %in% c("wu_natgen_CID4040", "wu_natgen_CID3838", "wu_natgen_CID4040", "BC393_Tumor", "BC392_Tumor") ~ "EFO:0009922"))
epi_rpca_new$disease_ontology_term_id <- "MONDO:0007254"

# All malignant?
epi_rpca_new$cell_type_ontology_term_id <- "CL:0001064"
epi_rpca_new$self_reported_ethnicity_ontology_term_id <- "unknown"
# Cell x gene will help you translate numerical age values to ontology ids
epi_rpca_new$development_stage_ontology_term_id <- epi_rpca_new$age
epi_rpca_new$sex_ontology_term_id <- "PATO:0000383"
epi_rpca_new$donor_id <- epi_rpca_new$donor
epi_rpca_new$suspension_type <- "cell"
epi_rpca_new$author_cell_type <- "Malignant"
epi_rpca_new@reductions$pca <- NULL
epi_rpca_new@reductions$umap.pca <- NULL
epi_rpca_new@reductions$integrated.rpca <- NULL
epi_rpca_new[["RNA"]] <- as(epi_rpca_new[["RNA"]], "Assay")
epi_rpca_new@meta.data <- epi_rpca_new@meta.data %>% dplyr::select(all_of(cellxgene_metadata_columns))
saveRDS(epi_rpca_new, file.path(PATH, "data/sc/cellxgene/epi.rds"))
sceasy::convertFormat(epi_rpca_new, from="seurat", to="anndata",
                      main_layer = "counts", transfer_layers="data",
                      drop_single_values = FALSE,
                      outFile=file.path(PATH, "data/sc/cellxgene/epi.h5ad"))
# Saving umap
combined_seurat_rpca <- readRDS("/restricted/projectnb/brcameta/brca_atlas/data/sc/combined_seurat_rpca.rds")
write.csv(combined_seurat_rpca@reductions$umap.rpca@cell.embeddings, file.path(PATH, "data/embeddings/all/rpca_umap.csv"), row.names = TRUE)

epi_rpca <- readRDS(file.path(PATH, "data/sc/cellxgene/epi.rds"))
strom_rpca <- readRDS(file.path(PATH, "data/sc/cellxgene/strom.rds"))
imm_rpca <- readRDS(file.path(PATH, "data/sc/cellxgene/imm.rds"))
combined_seurat <- scCustomize::Merge_Seurat_List(list(epi_rpca, imm_rpca, strom_rpca), project = "cellxgene")
saveRDS(combined_seurat, file.path(PATH, "data/sc/cellxgene/all.rds"))
