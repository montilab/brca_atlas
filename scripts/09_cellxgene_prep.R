library(Seurat)
library(tidyverse)
library(sceasy)
library("AnnotationDbi")
library("org.Hs.eg.db")
library(reticulate)
reticulate::use_condaenv("r-sceasy")
options(Seurat.object.assay.version = "v5")

PATH <- file.path(Sys.getenv("MLAB"), "projects/brcameta/brca_atlas/")

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
# # Add metadata information
# imm_rpca$organism_ontology_term_id <- "NCBITaxon:9606"
# imm_rpca$tissue_ontology_term_id <- "UBERON:0000310"
# imm_rpca$tissue_type <- "tissue"
# imm_rpca$assay_ontology_term_id <- with(imm_rpca@meta.data, 
#                                              case_when(batch == "wu_natgen_2021" ~ "EFO:0009899",
#                                                        batch == "wu_embo_2020" ~ "EFO:0009899",
#                                                        batch == "wu_genomemed_2021" ~ "EFO:0008995",
#                                                        batch == "pal_2021" ~ "EFO:0030003",
#                                                        batch == "qian_2020" ~ "EFO:0008995",
#                                                        batch == "gao_2021" ~ "EFO:0009899",
#                                                        batch == "tietscher_2023" ~ "EFO:0009922",
#                                                        batch == "liu_2023" ~ "EFO:0009922",
#                                                        batch == "wang_2024" ~ "EFO:0009899",
#                                                        batch == "bassez_2021" ~ "EFO:0030004"))
# imm_rpca$disease_ontology_term_id <- "MONDO:0007254"
# # Match Cell type Ontologies to annotations found in dataset
# imm_rpca$cell_type_ontology_term_id <- with(imm_rpca@meta.data,
#                                       case_when(`RNA_snn_res.0.4` == 0 ~ "CL:0000864",
#                                                 `RNA_snn_res.0.4` == 1 ~ "CL:0000623",
#                                                 `RNA_snn_res.0.4` == 2 ~ "CL:0000904",
#                                                 `RNA_snn_res.0.4` == 3 ~ "CL:0000907",
#                                                 `RNA_snn_res.0.4` == 4 ~ "CL:0011025",
#                                                 `RNA_snn_res.0.4` == 5 ~ "CL:0001048",
#                                                 `RNA_snn_res.0.4` == 6 ~ "CL:0000236",
#                                                 `RNA_snn_res.0.4` == 7 ~ "CL:0000990",
#                                                 `RNA_snn_res.0.4` == 8 ~ "CL:0000786",
#                                                 `RNA_snn_res.0.4` == 10 ~ "CL:0000235",
#                                                 `RNA_snn_res.0.4` == 11 ~ "CL:0000625",
#                                                 `RNA_snn_res.0.4` == 12 ~ "CL:4033069",
#                                                 `RNA_snn_res.0.4` == 13 ~ "CL:0000097",
#                                                 `RNA_snn_res.0.4` == 15 ~ "CL:0001058",
#                                                 `RNA_snn_res.0.4` == 17 ~ "CL:0001057",
#                                                 `RNA_snn_res.0.4` == 18 ~ "CL:4033076",
#                                                 `RNA_snn_res.0.4` == 19 ~ "CL:0000785",
#                                                 `RNA_snn_res.0.4` == 20 ~ "CL:0000623",
#                                                 `RNA_snn_res.0.4` == 21 ~ "CL:0002477"))
# imm_rpca$self_reported_ethnicity_ontology_term_id <- "unknown"
# # Cell x gene will help you translate numerical age values to ontology ids
# imm_rpca$development_stage_ontology_term_id <- imm_rpca$age
# imm_rpca$sex_ontology_term_id <- "PATO:0000383"
# imm_rpca$donor_id <- imm_rpca$donor
# imm_rpca$suspension_type <- "cell"
# imm_rpca$author_cell_type <- imm_rpca$cluster_annot
# imm_rpca@reductions$pca <- NULL
# imm_rpca@reductions$umap.pca <- NULL
# imm_rpca@reductions$integrated.rpca <- NULL
# imm_rpca[["RNA"]] <- as(imm_rpca[["RNA"]], "Assay")
# imm_rpca@meta.data <- imm_rpca@meta.data %>% dplyr::select(all_of(cellxgene_metadata_columns))
# sceasy::convertFormat(imm_rpca, from="seurat", to="anndata",
#                       main_layer = "counts", transfer_layers="data",
#                       drop_single_values = FALSE,
#                       outFile=file.path(PATH, "data/sc/imm_cellxgene.h5ad"))      
# 
# # 2. Strom Subset
# strom_rpca <- readRDS(file.path(PATH, "data/sc/strom_rpca_clean.rds"))
# 
# # Add metadata information
# strom_rpca$organism_ontology_term_id <- "NCBITaxon:9606"
# strom_rpca$tissue_ontology_term_id <- "UBERON:0000310"
# strom_rpca$tissue_type <- "tissue"
# strom_rpca$assay_ontology_term_id <- with(strom_rpca@meta.data, 
#                                         case_when(batch == "wu_natgen_2021" ~ "EFO:0009899",
#                                                   batch == "wu_embo_2020" ~ "EFO:0009899",
#                                                   batch == "wu_genomemed_2021" ~ "EFO:0008995",
#                                                   batch == "pal_2021" ~ "EFO:0030003",
#                                                   batch == "qian_2020" ~ "EFO:0008995",
#                                                   batch == "gao_2021" ~ "EFO:0009899",
#                                                   batch == "tietscher_2023" ~ "EFO:0009922",
#                                                   batch == "liu_2023" ~ "EFO:0009922",
#                                                   batch == "wang_2024" ~ "EFO:0009899",
#                                                   batch == "bassez_2021" ~ "EFO:0030004"))
# strom_rpca$disease_ontology_term_id <- "MONDO:0007254"
# # Match Cell type Ontologies to annotations found in dataset
# strom_rpca$cell_type_ontology_term_id <- with(strom_rpca@meta.data,
#                                             case_when(`RNA_snn_res.0.4` == 0 ~ "CL:0000057",
#                                                       `RNA_snn_res.0.4` == 1 ~ "CL:0002543",
#                                                       `RNA_snn_res.0.4` == 2 ~ "CL:0000057",
#                                                       `RNA_snn_res.0.4` == 3 ~ "CL:0000359",
#                                                       `RNA_snn_res.0.4` == 4 ~ "CL:0002144",
#                                                       `RNA_snn_res.0.4` == 5 ~ "CL:0000057",
#                                                       `RNA_snn_res.0.4` == 6 ~ "CL:0000669",
#                                                       `RNA_snn_res.0.4` == 7 ~ "CL:0000115",
#                                                       `RNA_snn_res.0.4` == 9 ~ "CL:1000413",
#                                                       `RNA_snn_res.0.4` == 10 ~ "CL:0000057",
#                                                       `RNA_snn_res.0.4` == 11 ~ "CL:0000669",
#                                                       `RNA_snn_res.0.4` == 12 ~ "CL:0002144",
#                                                       `RNA_snn_res.0.4` == 13 ~ "CL:0000115",
#                                                       `RNA_snn_res.0.4` == 14 ~ "CL:0002138",
#                                                       `RNA_snn_res.0.4` == 14 ~ "CL:0000115"))
# strom_rpca$self_reported_ethnicity_ontology_term_id <- "unknown"
# # Cell x gene will help you translate numerical age values to ontology ids
# strom_rpca$development_stage_ontology_term_id <- strom_rpca$age
# strom_rpca$sex_ontology_term_id <- "PATO:0000383"
# strom_rpca$donor_id <- strom_rpca$donor
# strom_rpca$suspension_type <- "cell"
# strom_rpca$author_cell_type <- strom_rpca$cluster_annot
# strom_rpca@reductions$pca <- NULL
# strom_rpca@reductions$umap.pca <- NULL
# strom_rpca@reductions$integrated.rpca <- NULL
# strom_rpca[["RNA"]] <- as(strom_rpca[["RNA"]], "Assay")
# strom_rpca@meta.data <- strom_rpca@meta.data %>% dplyr::select(all_of(cellxgene_metadata_columns))
# sceasy::convertFormat(strom_rpca, from="seurat", to="anndata",
#                       main_layer = "counts", transfer_layers="data",
#                       drop_single_values = FALSE,
#                       outFile=file.path(PATH, "data/sc/strom_cellxgene.h5ad")) 

# 3. Epi Subset
epi_rpca <- readRDS(file.path(PATH, "data/sc/epi_rpca_subset.rds"))
# Add metadata information
epi_rpca$organism_ontology_term_id <- "NCBITaxon:9606"
epi_rpca$tissue_ontology_term_id <- "UBERON:0000310"
epi_rpca$tissue_type <- "tissue"
epi_rpca$assay_ontology_term_id <- with(epi_rpca@meta.data, 
                                        case_when(batch == "wu_natgen_2021" ~ "EFO:0009899",
                                                  batch == "wu_embo_2020" ~ "EFO:0009899",
                                                  batch == "wu_genomemed_2021" ~ "EFO:0008995",
                                                  batch == "pal_2021" ~ "EFO:0030003",
                                                  batch == "qian_2020" ~ "EFO:0008995",
                                                  batch == "gao_2021" ~ "EFO:0009899",
                                                  batch == "tietscher_2023" ~ "EFO:0009922",
                                                  batch == "liu_2023" ~ "EFO:0009922",
                                                  batch == "wang_2024" ~ "EFO:0009899",
                                                  batch == "bassez_2021" ~ "EFO:0030004"))
epi_rpca$disease_ontology_term_id <- "MONDO:0007254"

# All malignant?
epi_rpca$cell_type_ontology_term_id <- "CL:0001064"
epi_rpca$self_reported_ethnicity_ontology_term_id <- "unknown"
# Cell x gene will help you translate numerical age values to ontology ids
epi_rpca$development_stage_ontology_term_id <- epi_rpca$age
epi_rpca$sex_ontology_term_id <- "PATO:0000383"
epi_rpca$donor_id <- epi_rpca$donor
epi_rpca$suspension_type <- "cell"
epi_rpca$author_cell_type <- "Malignant"
epi_rpca@reductions$pca <- NULL
epi_rpca@reductions$umap.pca <- NULL
epi_rpca@reductions$integrated.rpca <- NULL
epi_rpca[["RNA"]] <- as(epi_rpca[["RNA"]], "Assay")
epi_rpca@meta.data <- epi_rpca@meta.data %>% dplyr::select(all_of(cellxgene_metadata_columns))
sceasy::convertFormat(epi_rpca, from="seurat", to="anndata",
                      main_layer = "counts", transfer_layers="data",
                      drop_single_values = FALSE,
                      outFile=file.path(PATH, "data/sc/epi_cellxgene.h5ad")) 
# Saving umap
combined_seurat_rpca <- readRDS("/restricted/projectnb/brcameta/brca_atlas/data/sc/combined_seurat_rpca.rds")
write.csv(combined_seurat_rpca@reductions$umap.rpca@cell.embeddings, file.path(PATH, "data/embeddings/all/rpca_umap.csv"), row.names = TRUE)

