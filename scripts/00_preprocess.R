library(tidyverse)
library(Seurat)
library(Biobase)
library(BiocParallel)
library(SingleCellExperiment)
library(scDblFinder)
library("AnnotationDbi")
library("org.Hs.eg.db")
options(Seurat.object.assay.version = "v5")

DATA_PATH <- file.path(Sys.getenv("CBM"), "otherStudies/scRNAseq/brca")
PATH <- file.path(Sys.getenv("MLAB"), "projects/brcameta/brca_atlas")

# MAD-Filter Helper
mad_filter <- function(seurat_object, metric, n_mads) {
  metrics <- seurat_object@meta.data[[metric]]
  outlier <- (metrics < median(metrics) - n_mads*mad(metrics)) |
    (metrics > median(metrics) + n_mads*mad(metrics))
  if(all(is.na(outlier))) {
    outlier <- rep(FALSE, length.out=length(metrics))
  }
  return(outlier)
}

# 1. Wu Natgen 2021
wu_natgen_2021_exp <- Matrix::readMM(file.path(DATA_PATH, "wu_natgen_2021/raw/matrix.mtx"))
wu_natgen_2021_feats <- read.csv(file.path(DATA_PATH, "wu_natgen_2021/raw/features.tsv"), sep="\t", header=FALSE)
wu_natgen_2021_metadata <- read.csv(file.path(DATA_PATH, "wu_natgen_2021/raw/Whole_miniatlas_meta.csv"), skip=1)

rownames(wu_natgen_2021_exp) <- wu_natgen_2021_feats$V1
rownames(wu_natgen_2021_metadata) <- wu_natgen_2021_metadata$TYPE
rownames(wu_natgen_2021_feats) <- wu_natgen_2021_feats$V1
colnames(wu_natgen_2021_exp) <- rownames(wu_natgen_2021_metadata)

wu_natgen_2021_exp <- as(wu_natgen_2021_exp, "dgCMatrix")

# QC + Normalization
# All wu profiles already have been called with emptyDrops
wu_natgen_2021_seurat <- CreateSeuratObject(counts=wu_natgen_2021_exp,
                                            meta.data=wu_natgen_2021_metadata)
wu_natgen_2021_seurat$mitoRatio <- PercentageFeatureSet(object = wu_natgen_2021_seurat, pattern = "^MT-")
wu_natgen_2021_seurat$mitoRatio <- wu_natgen_2021_seurat@meta.data$mitoRatio / 100
wu_natgen_2021_seurat$nFeature_RNA %>% fivenum %>% unname
wu_natgen_2021_seurat$nCount_RNA %>% fivenum %>% unname
wu_natgen_2021_seurat$mitoRatio %>% fivenum %>% unname
# Cells have already been filtered. The paper describes the filters below.
# This subsetting doesn't remove any cells.
wu_natgen_2021_seurat <- subset(x = wu_natgen_2021_seurat,
                                subset= (nCount_RNA >= 250) &
                                  (nFeature_RNA >= 200) &
                                  (mitoRatio < 0.20))
wu_natgen_2021_seurat <- Seurat::NormalizeData(wu_natgen_2021_seurat)
wu_natgen_2021_sce <- Seurat::as.SingleCellExperiment(wu_natgen_2021_seurat)
wu_natgen_2021_sce <- scDblFinder::scDblFinder(wu_natgen_2021_sce,
                                               samples = "group",
                                               BPPARAM = MulticoreParam(3))

new_coldata <- S4Vectors::DataFrame(donor = paste0("wu_natgen_", wu_natgen_2021_sce$group),
                                    celltype_major = wu_natgen_2021_sce$group.1,
                                    celltype_minor = wu_natgen_2021_sce$group.3,
                                    subtype = wu_natgen_2021_sce$group.4,
                                    scDblFinder.score = wu_natgen_2021_sce$scDblFinder.score,
                                    scDblFinder.class = wu_natgen_2021_sce$scDblFinder.class)
rownames(new_coldata) <- colnames(wu_natgen_2021_sce)
colData(wu_natgen_2021_sce) <- new_coldata
saveRDS(wu_natgen_2021_sce, file=file.path(DATA_PATH, "wu_natgen_2021/processed/processed_sce.rds"))

# 2. Wu EMBO 2020
wu_embo_2020_exp_raw <- Matrix::readMM(file.path(DATA_PATH, "wu_embo_2020/raw/Wu_EMBOJ_count_matrix_sparse.mtx"))
wu_embo_2020_exp_genes <- read.delim(file.path(DATA_PATH, "wu_embo_2020/raw/Wu_EMBOJ_genes.tsv"), header=FALSE)
wu_embo_2020_exp_raw <- as(as.matrix(wu_embo_2020_exp_raw), "dgCMatrix")
wu_embo_2020_metadata <- read.csv(file.path(DATA_PATH, "wu_embo_2020/SCP1106/metadata/Wu_EMBO_metadata.csv"))
wu_embo_2020_metadata <- wu_embo_2020_metadata[-1,]

rownames(wu_embo_2020_metadata) <- wu_embo_2020_metadata$NAME
colnames(wu_embo_2020_exp_raw) <- wu_embo_2020_metadata$NAME
rownames(wu_embo_2020_exp_raw) <- wu_embo_2020_exp_genes$V1

# Normalizing using Seurat::LogNorm
wu_embo_2020_seurat <- CreateSeuratObject(counts=wu_embo_2020_exp_raw,
                                            meta.data=wu_embo_2020_metadata)

# Seems like raw data has already been QC'd (nUMI > 250, nGene > 200)
wu_embo_2020_seurat$mitoRatio <- PercentageFeatureSet(object = wu_embo_2020_seurat, pattern = "^MT-")
wu_embo_2020_seurat$mitoRatio <- wu_embo_2020_seurat@meta.data$mitoRatio / 100
wu_embo_2020_seurat$nFeature_RNA %>% fivenum
wu_embo_2020_seurat$nCount_RNA %>% fivenum
wu_embo_2020_seurat$mitoRatio %>% fivenum
# No cells were dropped after this
wu_embo_2020_seurat <- subset(x = wu_embo_2020_seurat,
                                subset= (nCount_RNA >= 250) &
                                  (nFeature_RNA >= 200) &
                                  (mitoRatio < 0.1))

wu_embo_2020_seurat <- Seurat::NormalizeData(wu_embo_2020_seurat)
wu_embo_2020_sce <- Seurat::as.SingleCellExperiment(wu_embo_2020_seurat)
wu_embo_2020_sce <- scDblFinder::scDblFinder(wu_embo_2020_sce,
                                               samples = "orig.ident",
                                               BPPARAM = MulticoreParam(3))

# Filtering coldata to only keep celltype identity
new_coldata <- DataFrame(donor = paste0("wu_embo_", wu_embo_2020_sce$orig.ident),
                         celltype_major = wu_embo_2020_sce$celltype_final,
                         celltype_minor = rep(NA, ncol(wu_embo_2020_sce)),
                         subtype = "TNBC",
                         scDblFinder.score = wu_embo_2020_sce$scDblFinder.score,
                         scDblFinder.class = wu_embo_2020_sce$scDblFinder.class)
rownames(new_coldata) <- colnames(wu_embo_2020_sce)
# This patient is treated
wu_embo_2020_sce <- wu_embo_2020_sce[,wu_embo_2020_sce$donor != "wu_embo_CID4513"]
colData(wu_embo_2020_sce) <- new_coldata
saveRDS(wu_embo_2020_sce, file=file.path(DATA_PATH, "wu_embo_2020/processed/processed_sce.rds"))

# 3. Wu Genomemed 2021
wu_genomemed_2021_exp_raw <- Matrix::readMM(file.path(DATA_PATH, "wu_genomemed_2021/Wu_etal_2021_allcells_raw_counts.mtx"))
wu_genomemed_2021_feats <- read.csv(file.path(DATA_PATH, "wu_genomemed_2021/Wu_etal_2021_allcells_genes.tsv"), header=FALSE)
wu_genomemed_2021_metadata <- read.csv(file.path(DATA_PATH, "wu_genomemed_2021/Wu_etal_2021_metadata.txt"), sep="\t", header=TRUE)
wu_genomemed_2021_metadata <- wu_genomemed_2021_metadata[-1,]

rownames(wu_genomemed_2021_metadata) <- wu_genomemed_2021_metadata$NAME
rownames(wu_genomemed_2021_feats) <- wu_genomemed_2021_feats$V1
rownames(wu_genomemed_2021_exp_raw) <- rownames(wu_genomemed_2021_exp_raw) <- wu_genomemed_2021_feats$V1
colnames(wu_genomemed_2021_exp_raw) <- colnames(wu_genomemed_2021_exp_raw) <- rownames(wu_genomemed_2021_metadata)
wu_genomemed_2021_exp_raw <- as(as.matrix(wu_genomemed_2021_exp_raw), "dgCMatrix")

wu_genomemed_2021_seurat <- CreateSeuratObject(counts=wu_genomemed_2021_exp_raw,
                                          meta.data=wu_genomemed_2021_metadata)
# QC + Normalization
# Same nUMI > 250, nGene > 200 filter already applied
wu_genomemed_2021_seurat$nUMI <- colSums(wu_genomemed_2021_exp_raw)
wu_genomemed_2021_seurat$nGene <- colSums(wu_genomemed_2021_exp_raw != 0)
wu_genomemed_2021_seurat$mitoRatio <- PercentageFeatureSet(object = wu_genomemed_2021_seurat, pattern = "^MT-")
wu_genomemed_2021_seurat$mitoRatio <- wu_genomemed_2021_seurat@meta.data$mitoRatio / 100
wu_genomemed_2021_seurat$nUMI %>% fivenum
wu_genomemed_2021_seurat$nGene %>% fivenum
wu_genomemed_2021_seurat$mitoRatio %>% fivenum
wu_genomemed_2021_seurat <- subset(x = wu_genomemed_2021_seurat,
                                   subset= (nCount_RNA >= 250) &
                                     (nFeature_RNA >= 200) &
                                     (mitoRatio < 0.20))

wu_genomemed_2021_seurat <- Seurat::NormalizeData(wu_genomemed_2021_seurat)
wu_genomemed_2021_sce <- Seurat::as.SingleCellExperiment(wu_genomemed_2021_seurat)
# breast cancer filter
breast_cancer_filter <- wu_genomemed_2021_sce$disease__ontology_label == "breast cancer"
wu_genomemed_2021_sce <- wu_genomemed_2021_sce[,breast_cancer_filter]
mouse_filter <- wu_genomemed_2021_sce$CellType != "Mouse 3T3 spike-in"
wu_genomemed_2021_sce <- wu_genomemed_2021_sce[,mouse_filter]
wu_genomemed_2021_sce <- scDblFinder::scDblFinder(wu_genomemed_2021_sce,
                                                  samples = "biosample_id")
# Filtering to matchable ensemblIDs
hgnc_symbols <- wu_genomemed_2021_sce %>% rownames
hgnc_symbols_new <- limma::alias2SymbolTable(hgnc_symbols, species="Hs")
old_new <- data.frame(old = hgnc_symbols, new = hgnc_symbols_new)
old_new$new <- with(old_new, if_else(is.na(new), true=old, false=new))
ensemblid <- mapIds(org.Hs.eg.db,
                    keys=old_new$new,
                    column="ENSEMBL",
                    keytype="SYMBOL",
                    multiVals="first")
ensemblid <- ensemblid[!is.na(ensemblid)]
ensembl_hgnc_df <- data.frame(ensembl = unname(ensemblid), hgnc=names(ensemblid))

rownames_filter <- old_new$new %in% ensembl_hgnc_df$hgnc
wu_genomemed_2021_sce <- wu_genomemed_2021_sce[rownames_filter,]

rowData(wu_genomemed_2021_sce) <- ensembl_hgnc_df
rownames(wu_genomemed_2021_sce) <- ensembl_hgnc_df$ensembl

# doublets & mouse filter
doublets_filter <- wu_genomemed_2021_sce$CellType != "Doublets"

# Filtering coldata to only keep celltype identity
new_coldata <- DataFrame(donor = paste0("wu_genomemed_", wu_genomemed_2021_sce$donor_id),
                         celltype_major = wu_genomemed_2021_sce$CellType,
                         celltype_minor = rep(NA, ncol(wu_genomemed_2021_sce)),
                         subtype = rep(NA, ncol(wu_genomemed_2021_sce)),
                         scDblFinder.score = wu_genomemed_2021_sce$scDblFinder.score,
                         scDblFinder.class = wu_genomemed_2021_sce$scDblFinder.class)
rownames(new_coldata) <- colnames(wu_genomemed_2021_sce)
colData(wu_genomemed_2021_sce) <- new_coldata
saveRDS(wu_genomemed_2021_sce, file=file.path(DATA_PATH, "wu_genomemed_2021/processed_sce.rds"))

# 4. Pal 2021
pal_paths <- Sys.glob(file.path(DATA_PATH, "pal_2021/GSM*.mtx"))
pal_barcodes_paths <- Sys.glob(file.path(DATA_PATH, "pal_2021/GSM*.tsv"))
pal_feats_path <- file.path(DATA_PATH, "pal_2021/GSE161529_features.tsv")
pal_sample_info <- read.delim(file.path(DATA_PATH, "pal_2021/GSE161529_sample_data.csv"), header=FALSE)
#
for(i in seq_along(pal_paths)) {
  pal_mtx <- Seurat::ReadMtx(pal_paths[[i]],
                             features = pal_feats_path,
                             cells = pal_barcodes_paths[[i]],
                             feature.column = 2 )
  pal_barcodes <- read.delim(pal_barcodes_paths[[i]], header=FALSE) %>% pull
  n_samples <- dim(pal_mtx)[[2]]
  pal_sample <- data.frame(sample = 1:n_samples, pheno=pal_sample_info[i,2])

  pal_seurat <- CreateSeuratObject(counts=pal_mtx,
                                   meta.data=pal_sample)
  pal_seurat$nUMI <- colSums(pal_mtx)
  pal_seurat$nGene <- colSums(pal_mtx != 0)
  pal_seurat$mitoRatio <- PercentageFeatureSet(object = pal_seurat, pattern = "^MT-")
  pal_seurat$mitoRatio <- pal_seurat@meta.data$mitoRatio / 100
  pal_seurat$nUMI %>% fivenum
  pal_seurat$nGene %>% fivenum
  pal_seurat$mitoRatio %>% fivenum
  pal_seurat$log1p_counts <- log1p(pal_seurat$nCount_RNA)
  pal_seurat$log1p_n_gene <- log1p(pal_seurat$nFeature_RNA)
  # At least 500 detected genes were generally required for each cell,
  # although the lower limit was reduced to 400 or 300 for some libraries.
  # No more than 20% mitochondrial reads were generally allowed per cell,
  # although the upper limit was increased as high as 40% for a small number of libraries.
  # pal_seurat <- subset(x = pal_seurat,
  #                                    subset= (nCount_RNA >= 500) &
  #                                      (nFeature_RNA >= 500) &
  #                                      (mitoRatio < 0.20))
  # Deciding to go with own filtering
  pal_seurat$low_quality <- mad_filter(pal_seurat, "log1p_n_gene", n_mads = 3) | mad_filter(pal_seurat, "log1p_counts", n_mads = 3) | mad_filter(pal_seurat, "mitoRatio", 3)
  cells_before <- dim(pal_seurat)[[2]]
  
  # # According to their methods section
  # seurat_obj <- subset(x = seurat_obj,
  #                      subset= (nUMI <= 75000) &
  #                        (nGene >= 200) &
  #                        (nGene < 7500) &
  #                        (mitoRatio < 0.20))
  
  # MAD Filtering
  pal_seurat <- subset(x = pal_seurat, subset = (low_quality != TRUE))
  
  cells_after <- dim(pal_seurat)[[2]]
  print(paste0("Removed ", cells_before - cells_after, " cells."))
  pal_seurat <- Seurat::NormalizeData(pal_seurat)

  pal_sce <- Seurat::as.SingleCellExperiment(pal_seurat)
  pal_sce <- scDblFinder::scDblFinder(pal_sce)

  saveRDS(pal_sce, file=file.path(DATA_PATH, paste0("pal_2021/processed/", pal_sample_info[i,2], ".rds")))
}

pal_paths <- Sys.glob(file.path(DATA_PATH, "pal_2021/processed/*Patient*.rds"))
pal_rds <- lapply(pal_paths, readRDS)
for (i in seq_along(pal_rds)) {
  rowData(pal_rds[[i]])$scDblFinder.selected <- NULL
}
pal_combined <- purrr::reduce(pal_rds, SingleCellExperiment::cbind)

metadata <- pal_combined$pheno
tumor_filter <- grepl("tumour.*", metadata)
lymph_filter <- !grepl("Lymph.*", metadata)
tumor_and_not_lymph_filter <- tumor_filter & lymph_filter
batch <- str_extract(metadata, "Patient.+")
subtype <- str_extract(metadata, "(.+)tumour", group=1)
new_coldata <- DataFrame(donor = paste0("pal_", batch),
                         celltype_major = rep(NA, ncol(pal_combined)),
                         celltype_minor = rep(NA, ncol(pal_combined)),
                         subtype = subtype,
                         batch = batch,
                         scDblFinder.score = pal_combined$scDblFinder.score,
                         scDblFinder.class = pal_combined$scDblFinder.class)
rownames(new_coldata) <- colnames(pal_combined)
colData(pal_combined) <- new_coldata

saveRDS(pal_combined, file=file.path(DATA_PATH, "pal_2021/processed/combined.rds"))

pal_cancer <- pal_combined[,tumor_and_not_lymph_filter]
saveRDS(pal_cancer, file=file.path(DATA_PATH, "pal_2021/processed/combined_cancer.rds"))

# 5. Qian 2020 (Raw version)
qian_counts <- Matrix::readMM(file.path(DATA_PATH, "qian_2020/export/BC_counts/matrix.mtx"))
qian_barcodes <- read.delim(file.path(DATA_PATH, "qian_2020/export/BC_counts/barcodes.tsv"), header = FALSE)
qian_genes <- read.delim(file.path(DATA_PATH, "qian_2020/export/BC_counts/genes.tsv"), header = FALSE)
qian_metadata <- read.csv(file.path(DATA_PATH, "qian_2020/2103-Breastcancer_metadata.csv.gz"))
stopifnot(all(qian_barcodes == qian_metadata$Cell)) #Barcode information is same as metadata
stopifnot(all(qian_genes$V1 == qian_genes$V2))
rownames(qian_metadata) <- qian_metadata$Cell
rownames(qian_counts) <- qian_genes$V1
colnames(qian_counts) <- rownames(qian_metadata)
qian_counts <- as(as.matrix(qian_counts), "dgCMatrix")
# Normalizing using Seurat::LogNorm
qian_seurat <- CreateSeuratObject(counts=qian_counts,
                                  meta.data=qian_metadata)
# QC + Normalization
# Same as paper
qian_seurat$nUMI <- colSums(qian_counts)
qian_seurat$nGene <- colSums(qian_counts != 0)
all(qian_seurat$nUMI == qian_seurat$nCount_RNA)
all(qian_seurat$nGene == qian_seurat$nFeature_RNA)
qian_seurat$mitoRatio <- PercentageFeatureSet(object = qian_seurat, pattern = "^MT-")
qian_seurat$mitoRatio <- qian_seurat@meta.data$mitoRatio / 100
qian_seurat$nUMI %>% fivenum
qian_seurat$nGene %>% fivenum
qian_seurat$mitoRatio %>% fivenum
qian_seurat <- subset(x = qian_seurat,
                      subset= (nCount_RNA >= 401) &
                        (nFeature_RNA >= 201) &
                        (nFeature_RNA <= 6000) &
                        (mitoRatio < 0.25))
qian_seurat <- Seurat::NormalizeData(qian_seurat)
qian_sce <- Seurat::as.SingleCellExperiment(qian_seurat)
qian_sce <- scDblFinder::scDblFinder(qian_sce, samples = 'orig.ident')
# Filtering coldata to only keep celltype identity
new_coldata <- DataFrame(donor = paste0("qian_", qian_sce$PatientNumber),
                         celltype_major = qian_sce$CellType,
                         celltype_minor = rep(NA, ncol(qian_sce)),
                         subtype = rep(NA, ncol(qian_sce)),
                         scDblFinder.score = qian_sce$scDblFinder.score,
                         scDblFinder.class = qian_sce$scDblFinder.class)
new_coldata$celltype_major[new_coldata$celltype_major == ""] <- NA
rownames(new_coldata) <- colnames(qian_sce)

colData(qian_sce) <- new_coldata
saveRDS(qian_sce, file=file.path(DATA_PATH, "qian_2020/processed_sce_v2.rds"))
# 5. Qian 2020 (3CE version)
# qian_counts <- Matrix::readMM(file.path(DATA_PATH, "qian_2020/raw/Exp_data_UMIcounts.mtx"))
# qian_metadata <- read.csv(file.path(DATA_PATH, "qian_2020/raw/Cells.csv"))
# qian_feats <- read.delim(file.path(DATA_PATH, "qian_2020/raw/Genes.txt"), header=FALSE)
# 
# rownames(qian_metadata) <- qian_metadata$cell_name
# rownames(qian_counts) <- qian_feats$V1
# colnames(qian_counts) <- rownames(qian_metadata)
# 
# qian_counts <- as(as.matrix(qian_counts), "dgCMatrix")
# 
# # Normalizing using Seurat::LogNorm
# qian_seurat <- CreateSeuratObject(counts=qian_counts,
#                                   meta.data=qian_metadata)
# # QC + Normalization
# # Same as paper
# qian_seurat$nUMI <- colSums(qian_counts)
# qian_seurat$nGene <- colSums(qian_counts != 0)
# qian_seurat$mitoRatio <- PercentageFeatureSet(object = qian_seurat, pattern = "^MT-")
# qian_seurat$mitoRatio <- qian_seurat@meta.data$mitoRatio / 100
# qian_seurat$nUMI %>% fivenum
# qian_seurat$nGene %>% fivenum
# qian_seurat$mitoRatio %>% fivenum
# qian_seurat <- subset(x = qian_seurat,
#                       subset= (nCount_RNA >= 401) &
#                         (nFeature_RNA >= 201) &
#                         (nFeature_RNA <= 6000) &
#                         (mitoRatio < 0.25))
# qian_seurat <- Seurat::NormalizeData(qian_seurat)
# qian_sce <- Seurat::as.SingleCellExperiment(qian_seurat)
# qian_sce <- scDblFinder::scDblFinder(qian_sce, samples = 'orig.ident')
# #
# # Filtering to matchable ensemblIDs
# # hgnc_symbols <- qian_sce %>% rownames
# # hgnc_symbols_new <- limma::alias2SymbolTable(hgnc_symbols, species="Hs")
# # old_new <- data.frame(old = hgnc_symbols, new = hgnc_symbols_new)
# # old_new$new <- with(old_new, if_else(is.na(new), true=old, false=new))
# # ensemblid <- mapIds(org.Hs.eg.db,
# #                     keys=old_new$new,
# #                     column="ENSEMBL",
# #                     keytype="SYMBOL",
# #                     multiVals="first")
# # ensemblid <- ensemblid[!is.na(ensemblid)]
# # ensembl_hgnc_df <- data.frame(ensembl = unname(ensemblid), hgnc=names(ensemblid))
# #
# # rownames_filter <- old_new$new %in% ensembl_hgnc_df$hgnc
# # qian_sce <- qian_sce[rownames_filter,]
# #
# # rowData(qian_sce) <- ensembl_hgnc_df
# # rownames(qian_sce) <- ensembl_hgnc_df$ensembl
# #
# # Filtering coldata to only keep celltype identity
# new_coldata <- DataFrame(donor = paste0("qian_", qian_sce$cell_name %>% str_extract("(.+)_", group=1)),
#                          celltype_major = qian_sce$cell_type,
#                          celltype_minor = rep(NA, ncol(qian_sce)),
#                          subtype = rep(NA, ncol(qian_sce)),
#                          scDblFinder.score = qian_sce$scDblFinder.score,
#                          scDblFinder.class = qian_sce$scDblFinder.class)
# new_coldata$celltype_major[new_coldata$celltype_major == ""] <- NA
# rownames(new_coldata) <- colnames(qian_sce)
# 
# colData(qian_sce) <- new_coldata
# saveRDS(qian_sce, file=file.path(DATA_PATH, "qian_2020/processed_sce.rds"))

# 6. Gao 2021
gao_exp <- Matrix::readMM(file.path(DATA_PATH, "gao_2021/Breast/Exp_data_UMIcounts.mtx"))
gao_features <- read.delim(file.path(DATA_PATH, "gao_2021/Breast/Genes.txt"), header=FALSE)
# Though there's 5 TNBC samples, only 3 are included in the cells metadata below
gao_metadata <- read.csv(file.path(DATA_PATH, "gao_2021/Meta-data.csv"))
gao_cells <- read.csv(file.path(DATA_PATH, "gao_2021/Breast/Cells.csv"))

gao_cells <- merge(gao_metadata, gao_cells, by = "sample") %>%
  dplyr::select(cell_name, patient, cancer_type, sample, cell_type, histology)
gao_cells <- column_to_rownames(gao_cells, var="cell_name")

rownames(gao_exp) <- gao_features$V1
colnames(gao_exp) <- rownames(gao_cells)
gao_exp <- as(as.matrix(gao_exp), "dgCMatrix")

gao_seurat <- CreateSeuratObject(counts=gao_exp,
                                 meta.data=gao_cells)
gao_seurat <- gao_seurat[,gao_cells$cancer_type == "Breast cancer"]
gao_seurat$nUMI <- colSums(gao_exp)
gao_seurat$nGene <- colSums(gao_exp != 0)
gao_seurat$mitoRatio <- PercentageFeatureSet(object = gao_seurat, pattern = "^MT-")
gao_seurat$mitoRatio <- gao_seurat@meta.data$mitoRatio / 100
gao_seurat$nUMI %>% fivenum
gao_seurat$nGene %>% fivenum
gao_seurat$mitoRatio %>% fivenum
#genes_percent_expression_filter <- (rowMeans(gao_exp>0)*100) > 5
gao_seurat <- subset(x = gao_seurat,
                     subset= (nCount_RNA >= 250) &
                     (nFeature_RNA >= 200))
#gao_seurat <- gao_seurat[genes_percent_expression_filter,]
gao_seurat <- Seurat::NormalizeData(gao_seurat)
gao_sce <- Seurat::as.SingleCellExperiment(gao_seurat)
gao_sce <- scDblFinder::scDblFinder(gao_sce, sample = "sample")

# Filtering to matchable ensemblIDs
hgnc_symbols <- gao_sce %>% rownames
hgnc_symbols_new <- limma::alias2SymbolTable(hgnc_symbols, species="Hs")
old_new <- data.frame(old = hgnc_symbols, new = hgnc_symbols_new)
old_new$new <- with(old_new, if_else(is.na(new), true=old, false=new))
ensemblid <- mapIds(org.Hs.eg.db,
                    keys=old_new$new,
                    column="ENSEMBL",
                    keytype="SYMBOL",
                    multiVals="first")
ensemblid <- ensemblid[!is.na(ensemblid)]
ensembl_hgnc_df <- data.frame(ensembl = unname(ensemblid), hgnc=names(ensemblid))

rownames_filter <- old_new$new %in% ensembl_hgnc_df$hgnc
gao_sce <- gao_sce[rownames_filter,]

rowData(gao_sce) <- ensembl_hgnc_df
rownames(gao_sce) <- ensembl_hgnc_df$ensembl

# Cell type labels
new_coldata <- DataFrame(donor = paste0("gao_", gao_sce$patient),
                         celltype_major = gao_sce$cell_type,
                         celltype_minor = rep(NA, ncol(gao_sce)),
                         subtype = gao_sce$histology,
                         scDblFinder.score = gao_sce$scDblFinder.score,
                         scDblFinder.class = gao_sce$scDblFinder.class)
rownames(new_coldata) <- colnames(gao_sce)
colData(gao_sce) <- new_coldata

gao_sce$celltype_major[gao_sce$celltype_major == ""] <- NA
saveRDS(gao_sce, file=file.path(DATA_PATH, "gao_2021/processed_sce.rds"))

# 7. Tietscher 2023
matrix_paths <- Sys.glob(file.path(DATA_PATH, "tietscher_2023/*_matrix.txt"))
metadata_paths <- Sys.glob(file.path(DATA_PATH, "tietscher_2023/*_metadata.txt"))
clinical_data <- read.csv(file.path(DATA_PATH, "tietscher_2023/tietscher_clinical.csv"))
for(i in seq_along(matrix_paths)) {
  # Confirm if metadata and count matrix are matching
  sample_name <- str_match(matrix_paths[i], "tietscher_2023/(.*?)_singlecell")[,2]
  sample_name_metadata <- str_match(metadata_paths[i], "tietscher_2023/(.*?)_complete")[,2]
  stopifnot(sample_name == sample_name_metadata)
  metadata <- read.csv(metadata_paths[i], sep="\t")
  rownames(metadata) <- metadata$cellID
  stopifnot(sample_name == unique(metadata$sample))

  # Find subtype
  sample_name_clinical <- str_replace(sample_name, "T", "")
  subtype <- clinical_data %>% dplyr::filter(`Patient.ID` == sample_name_clinical) %>% pull(Clinical.Subtype)

  # Load count matrix
  count_matrix <- read.csv(matrix_paths[i], sep="\t")
  # For some samples there's less cells in metadata than in count matrix
  # Perhaps some were filtered already
  if (ncol(count_matrix) > nrow(metadata)) {
    print(sample_name)
    print("More cells in count matrix than metadata.")
  }
  na_filter <- rownames(count_matrix)[!is.na(rowSums(count_matrix))]
  colnames(count_matrix) <- paste0(sample_name, "_", str_replace(colnames(count_matrix), ".1", ""))
  cells_filter <- colnames(count_matrix) %in% rownames(metadata)
  count_matrix <- count_matrix[na_filter,cells_filter]
  stopifnot(ncol(count_matrix) == nrow(metadata))

  if (metadata$sample[1] %in% c("TBB165", "TBB171")) {
    print("Skipping neoadjuvant samples")
    next
  }

  count_matrix <- as(as.matrix(count_matrix), "dgCMatrix")

  seurat_obj <- CreateSeuratObject(counts=count_matrix,
                                   meta.data=metadata)

  seurat_obj$mitoRatio <- PercentageFeatureSet(object = seurat_obj, pattern = "^MT-")
  seurat_obj$mitoRatio <- seurat_obj@meta.data$mitoRatio / 100
  print(seurat_obj$mitoRatio %>% fivenum)
  print(seurat_obj$nFeature_RNA %>% fivenum)
  print(seurat_obj$nCount_RNA %>% fivenum)
  seurat_obj$log1p_counts <- log1p(seurat_obj$nCount_RNA)
  seurat_obj$log1p_n_gene <- log1p(seurat_obj$nFeature_RNA)
  seurat_obj$celltype_major <- seurat_obj$cell_type
  seurat_obj$celltype_minor <- NA
  seurat_obj$donor <- seurat_obj$sample
  seurat_obj$subtype <- subtype
  # Seems like uploaded data has yet to be filtered
  seurat_obj$low_quality <- mad_filter(seurat_obj, "log1p_n_gene", n_mads = 3) | mad_filter(seurat_obj, "log1p_counts", n_mads = 3) | mad_filter(seurat_obj, "mitoRatio", 3)
  cells_before <- dim(seurat_obj)[[2]]

  # # According to their methods section
  # seurat_obj <- subset(x = seurat_obj,
  #                      subset= (nUMI <= 75000) &
  #                        (nGene >= 200) &
  #                        (nGene < 7500) &
  #                        (mitoRatio < 0.20))

  # MAD Filtering
  seurat_obj <- subset(x = seurat_obj, subset = (low_quality != TRUE))

  cells_after <- dim(seurat_obj)[[2]]
  print(paste0("Removed ", cells_before - cells_after, " cells."))
  seurat_obj <- Seurat::NormalizeData(seurat_obj)

  sce_obj <- Seurat::as.SingleCellExperiment(seurat_obj)
  sce_obj <- scDblFinder::scDblFinder(sce_obj)

  saveRDS(sce_obj, file=file.path(DATA_PATH, paste0("tietscher_2023/processed/", sample_name, "_sce.rds")))
}

processed_paths <- Sys.glob(file.path(DATA_PATH, "tietscher_2023/processed/*_sce.rds"))
tietscher_rds <- lapply(processed_paths, readRDS)
for (i in seq_along(tietscher_rds)) {
  rowData(tietscher_rds[[i]])$scDblFinder.selected <- NULL
}
tietscher_common_features <- purrr::reduce(lapply(tietscher_rds, rownames), intersect)
tietscher_rds <- lapply(tietscher_rds, function(x) x[rownames(x) %in% tietscher_common_features,])
tietscher_combined <- purrr::reduce(tietscher_rds, SingleCellExperiment::cbind)
tietscher_combined$batch <- "tietscher_2023"
saveRDS(tietscher_combined, file=file.path(DATA_PATH, "tietscher_2023/processed/combined.rds"))

# 8. Bassez 2021
bassez_cohort1_counts <- readRDS(file.path(DATA_PATH, "bassez_2021/1863-counts_cells_cohort1.rds"))
bassez_cohort1_md <- read.csv(file.path(DATA_PATH, "bassez_2021/1872-BIOKEY_metaData_cohort1_web.csv"))
seurat_obj <- CreateSeuratObject(counts=bassez_cohort1_counts,
                                 meta.data=bassez_cohort1_md)
seurat_obj$mitoRatio <- PercentageFeatureSet(object = seurat_obj, pattern = "^MT-")
seurat_obj$mitoRatio <- seurat_obj@meta.data$mitoRatio / 100
seurat_obj$nCount_RNA %>% fivenum
seurat_obj$nFeature_RNA %>% fivenum
seurat_obj$mitoRatio %>% fivenum
# # Bassez is already filtered.
# seurat_obj <- subset(x = seurat_obj,
#                      subset= (nCount_RNA > 400) &
#                        (nFeature_RNA >= 200) &
#                        (nFeature_RNA < 6000) &
#                        (mitoRatio < 0.15))
seurat_obj <- subset(x = seurat_obj,
                     subset = (timepoint == "Pre"))
seurat_obj <- Seurat::NormalizeData(seurat_obj)
seurat_obj$celltype_major <- seurat_obj$cellType
seurat_obj$celltype_minor <- NA
seurat_obj$donor <- seurat_obj$patient_id
seurat_obj$subtype <- seurat_obj$BC_type
sce_obj <- Seurat::as.SingleCellExperiment(seurat_obj)
sce_obj <- scDblFinder::scDblFinder(sce_obj, samples = "patient_id")

saveRDS(sce_obj, file=file.path(DATA_PATH, "bassez_2021/processed_cohort1.rds"))

# 9. Xu 2021
xu_paths <- Sys.glob(file.path(DATA_PATH, "xu_2021/GSM*.txt"))
xu_primary_paths <- xu_paths[c(1,4,7,10,13)]
xu_subtype <- c("ER+", "HER2+", "HER2+", "TNBC", "TNBC")
xu_age <- c(54, 56, 68, 44, 40)
for(i in seq_along(xu_primary_paths)) {
  sample_name <- str_match(xu_primary_paths[i], "xu_2021/(.*?)_")
  sample_name <- sample_name[,2]

  count_matrix <- read.csv(xu_primary_paths[i], sep="\t")

  count_matrix <- as(as.matrix(count_matrix), "dgCMatrix")

  seurat_obj <- CreateSeuratObject(counts=count_matrix)

  seurat_obj$mitoRatio <- PercentageFeatureSet(object = seurat_obj, pattern = "^MT-")
  seurat_obj$mitoRatio <- seurat_obj@meta.data$mitoRatio / 100
  print(seurat_obj$nCount_RNA %>% fivenum)
  print(seurat_obj$nFeature_RNA %>% fivenum)
  print(seurat_obj$mitoRatio %>% fivenum)
  seurat_obj$subtype <- xu_subtype[i]
  seurat_obj$donor <- sample_name
  seurat_obj$age <- xu_age[i]
  seurat_obj$log1p_counts <- log1p(seurat_obj$nCount_RNA)
  seurat_obj$log1p_n_gene <- log1p(seurat_obj$nFeature_RNA)

  seurat_obj$low_quality <- mad_filter(seurat_obj, "log1p_n_gene", n_mads = 3) | mad_filter(seurat_obj, "log1p_counts", n_mads = 3) | mad_filter(seurat_obj, "mitoRatio", 3)
  # # According to their methods section
  # seurat_obj <- subset(x = seurat_obj,
  #                      subset= (nCount_RNA > 200) &
  #                        (nFeature_RNA > 100) &
  #                        (mitoRatio < 0.5))
  seurat_obj <- subset(x = seurat_obj, subset = (low_quality != TRUE))
  seurat_obj <- Seurat::NormalizeData(seurat_obj)

  sce_obj <- Seurat::as.SingleCellExperiment(seurat_obj)
  sce_obj <- scDblFinder::scDblFinder(sce_obj)
  saveRDS(sce_obj, file=file.path(DATA_PATH, paste0("xu_2021/processed/", sample_name, "_sce.rds")))
}
processed_paths <- Sys.glob(file.path(DATA_PATH, "xu_2021/processed/*_sce.rds"))
xu_rds <- lapply(processed_paths, readRDS)
for (i in seq_along(xu_rds)) {
  rowData(xu_rds[[i]])$scDblFinder.selected <- NULL
}
xu_common_features <- purrr::reduce(lapply(xu_rds, rownames), intersect)
xu_rds <- lapply(xu_rds, function(x) x[rownames(x) %in% xu_common_features,])
xu_combined <- purrr::reduce(xu_rds, SingleCellExperiment::cbind)
xu_combined$batch <- "xu_2021"
xu_combined$celltype_major <- NA
xu_combined$celltype_minor <- NA
saveRDS(xu_combined, file=file.path(DATA_PATH, "xu_2021/processed/combined.rds"))

# 10. Barkley 2022
barkley_matrices <- Sys.glob(file.path(DATA_PATH, "barkley_2022/*.mtx.gz"))
barkley_genes <- Sys.glob(file.path(DATA_PATH, "barkley_2022/*.tsv.gz"))

for(i in seq_along(barkley_matrices)) {
  sample_name <- str_match(barkley_matrices[i], "barkley_2022/.+(BRCA.)_")
  sample_name <- sample_name[,2]

  count_matrix <- Matrix::readMM(barkley_matrices[i])
  genes <- read_tsv(barkley_genes[[i]])
  count_matrix <- as(as.matrix(count_matrix), "dgCMatrix")
  rownames(count_matrix) <- genes$name

  seurat_obj <- CreateSeuratObject(counts=count_matrix)

  seurat_obj$mitoRatio <- PercentageFeatureSet(object = seurat_obj, pattern = "^MT-")
  seurat_obj$mitoRatio <- seurat_obj@meta.data$mitoRatio / 100
  print(seurat_obj$nCount_RNA %>% fivenum)
  print(seurat_obj$nFeature_RNA %>% fivenum)
  print(seurat_obj$mitoRatio %>% fivenum)
  seurat_obj$subtype <- NA
  seurat_obj$donor <- sample_name
  seurat_obj$age <- NA
  seurat_obj$log1p_counts <- log1p(seurat_obj$nCount_RNA)
  seurat_obj$log1p_n_gene <- log1p(seurat_obj$nFeature_RNA)

  seurat_obj$low_quality <- mad_filter(seurat_obj, "log1p_n_gene", n_mads = 3) | mad_filter(seurat_obj, "log1p_counts", n_mads = 3) | mad_filter(seurat_obj, "mitoRatio", 3)
  # # According to their methods section
  # seurat_obj <- subset(x = seurat_obj,
  #                      subset= (nCount_RNA > 200) &
  #                        (nFeature_RNA > 100) &
  #                        (mitoRatio < 0.5))
  seurat_obj <- subset(x = seurat_obj, subset = (low_quality != TRUE))
  seurat_obj <- Seurat::NormalizeData(seurat_obj)

  sce_obj <- Seurat::as.SingleCellExperiment(seurat_obj)
  sce_obj <- scDblFinder::scDblFinder(sce_obj)
  saveRDS(sce_obj, file=file.path(DATA_PATH, paste0("barkley_2022/processed/", sample_name, "_sce.rds")))
}
processed_paths <- Sys.glob(file.path(DATA_PATH, "barkley_2022/processed/*_sce.rds"))
barkley_rds <- lapply(processed_paths, readRDS)
for (i in seq_along(barkley_rds)) {
  rowData(barkley_rds[[i]])$scDblFinder.selected <- NULL
}
barkley_common_features <- purrr::reduce(lapply(barkley_rds, rownames), intersect)
barkley_rds <- lapply(barkley_rds, function(x) x[rownames(x) %in% barkley_common_features,])
barkley_combined <- purrr::reduce(barkley_rds, SingleCellExperiment::cbind)
barkley_combined$batch <- "barkley_2022"
barkley_combined$celltype_major <- NA
barkley_combined$celltype_minor <- NA
saveRDS(barkley_combined, file=file.path(DATA_PATH, "barkley_2022/processed/combined.rds"))

# 11. Liu
liu_counts <- Matrix::readMM(file.path(DATA_PATH, "liu_2023/GSE225600_sc_matrix.mtx.gz"))
liu_barcodes <- read.delim(file.path(DATA_PATH, "liu_2023/GSE225600_sc_barcodes.tsv.gz"), header = FALSE)
liu_genes <- read.delim(file.path(DATA_PATH, "liu_2023/GSE225600_sc_features.tsv.gz"), header = FALSE)
#Metadata from supp3 doesn't have matching cell barcodes
#liu_metadata <- readxl::read_excel(file.path(DATA_PATH, "liu_2023/ADVS-10-2205395-s006.xlsx"), sheet = "anniotation")
liu_barcodes$cellid <- lapply(str_split(liu_barcodes$V1, "-"), function(x) x[[1]]) %>% unlist
liu_barcodes$tumorid <- lapply(str_split(liu_barcodes$V1, "-"), function(x) x[[2]]) %>% unlist
rownames(liu_barcodes) <- liu_barcodes$V1
liu_barcodes$V1 <- NULL
rownames(liu_counts) <- liu_genes$V1
colnames(liu_counts) <- rownames(liu_barcodes)
liu_counts <- liu_counts[!duplicated(liu_genes$V1),]
liu_counts <- as(as.matrix(liu_counts), "dgCMatrix")
# Normalizing using Seurat::LogNorm
liu_seurat <- CreateSeuratObject(counts=liu_counts,
                                 meta.data=liu_barcodes)
# QC + Normalization
# Same as paper
liu_seurat$nUMI <- colSums(liu_counts)
liu_seurat$nGene <- colSums(liu_counts != 0)
liu_seurat$mitoRatio <- PercentageFeatureSet(object = liu_seurat, pattern = "^MT-")
liu_seurat$mitoRatio <- liu_seurat@meta.data$mitoRatio / 100
liu_seurat$nUMI %>% fivenum
liu_seurat$nGene %>% fivenum
liu_seurat$mitoRatio %>% fivenum
liu_seurat$log1p_counts <- log1p(liu_seurat$nUMI)
liu_seurat$log1p_n_gene <- log1p(liu_seurat$nGene)
# Decided to go with own filtering since raw dataset was not filtered
liu_seurat$low_quality <- mad_filter(liu_seurat, "log1p_n_gene", n_mads = 3) | mad_filter(liu_seurat, "log1p_counts", n_mads = 3) | mad_filter(liu_seurat, "mitoRatio", 3)
cells_before <- dim(liu_seurat)[[2]]

# MAD Filtering
liu_seurat <- subset(x = liu_seurat, subset = (low_quality != TRUE))

cells_after <- dim(liu_seurat)[[2]]
print(paste0("Removed ", cells_before - cells_after, " cells."))
liu_seurat <- Seurat::NormalizeData(liu_seurat)
liu_sce <- Seurat::as.SingleCellExperiment(liu_seurat)
liu_sce <- scDblFinder::scDblFinder(liu_sce, samples = 'orig.ident')
liu_sce <- liu_sce[,str_detect(liu_sce$tumorid, "T")]
# Filtering coldata to only keep celltype identity
new_coldata <- DataFrame(donor = liu_sce$tumorid,
                         celltype_major = rep(NA, ncol(liu_sce)),
                         celltype_minor = rep(NA, ncol(liu_sce)),
                         subtype = rep(NA, ncol(liu_sce)),
                         scDblFinder.score = liu_sce$scDblFinder.score,
                         scDblFinder.class = liu_sce$scDblFinder.class)
rownames(new_coldata) <- colnames(liu_sce)

colData(liu_sce) <- new_coldata
saveRDS(liu_sce, file=file.path(DATA_PATH, "liu_2023/processed_sce.rds"))

# 12. Wang
wang_matrices <- Sys.glob(file.path(DATA_PATH, "wang_2024/ER/*.mtx.gz"))
wang_genes <- Sys.glob(file.path(DATA_PATH, "wang_2024/ER/*features.tsv.gz"))
wang_barcodes <- Sys.glob(file.path(DATA_PATH, "wang_2024/ER/*barcodes.tsv.gz"))

for(i in seq_along(wang_matrices)) {
  sample_name <- str_match(wang_matrices[i], "ER/.+(BC.+)_")
  sample_name <- sample_name[,2]
  
  count_matrix <- Matrix::readMM(wang_matrices[i])
  genes <- read.delim(wang_genes[[i]], header = FALSE)
  barcodes <- read.delim(wang_barcodes[[i]], header = FALSE)
  count_matrix <- as(as.matrix(count_matrix), "dgCMatrix")
  rownames(count_matrix) <- genes$V2
  colnames(count_matrix) <- barcodes$V1
  count_matrix <- count_matrix[!duplicated(genes$V2),]
  seurat_obj <- CreateSeuratObject(counts=count_matrix)
  
  seurat_obj$nUMI <- colSums(count_matrix)
  seurat_obj$nGene <- colSums(count_matrix != 0)
  seurat_obj$mitoRatio <- PercentageFeatureSet(object = seurat_obj, pattern = "^MT-")
  seurat_obj$mitoRatio <- seurat_obj@meta.data$mitoRatio / 100
  print(seurat_obj$nUMI %>% fivenum)
  print(seurat_obj$nGene %>% fivenum)
  print(seurat_obj$mitoRatio %>% fivenum)
  next
  seurat_obj$log1p_counts <- log1p(seurat_obj$nCount_RNA)
  seurat_obj$log1p_n_gene <- log1p(seurat_obj$nFeature_RNA)
  seurat_obj$donor <- sample_name
  # Seems like nUMI was already filtered
  seurat_obj$low_quality <- mad_filter(seurat_obj, "log1p_n_gene", n_mads = 3) | mad_filter(seurat_obj, "mitoRatio", 3)
  # # According to their methods section
  # seurat_obj <- subset(x = seurat_obj,
  #                      subset= (nGene > 200) &
  #                        (nGene < 2500)
  #                        (mitoRatio < 0.2))
  seurat_obj <- subset(x = seurat_obj, subset = (low_quality != TRUE))
  seurat_obj <- Seurat::NormalizeData(seurat_obj)
  
  sce_obj <- Seurat::as.SingleCellExperiment(seurat_obj)
  sce_obj <- scDblFinder::scDblFinder(sce_obj)
  saveRDS(sce_obj, file=file.path(DATA_PATH, paste0("wang_2024/processed/ER/", sample_name, "_sce.rds")))
}
processed_paths <- Sys.glob(file.path(DATA_PATH, "wang_2024/processed/ER/*_sce.rds"))
wang_rds <- lapply(processed_paths, readRDS)
for (i in seq_along(wang_rds)) {
  rowData(wang_rds[[i]])$scDblFinder.selected <- NULL
}
wang_common_features <- purrr::reduce(lapply(wang_rds, rownames), intersect)
wang_rds <- lapply(wang_rds, function(x) x[rownames(x) %in% wang_common_features,])
wang_combined <- purrr::reduce(wang_rds, SingleCellExperiment::cbind)
wang_combined$batch <- "wang_2024"
wang_combined$celltype_major <- NA
wang_combined$celltype_minor <- NA
saveRDS(wang_combined, file=file.path(DATA_PATH, "wang_2024/processed/ER/combined.rds"))

# 13. Wang 2
wang_matrices <- Sys.glob(file.path(DATA_PATH, "wang_2024/TNBC/*.mtx.gz"))
wang_genes <- Sys.glob(file.path(DATA_PATH, "wang_2024/TNBC/*features.tsv.gz"))
wang_barcodes <- Sys.glob(file.path(DATA_PATH, "wang_2024/TNBC/*barcodes.tsv.gz"))

for(i in seq_along(wang_matrices)) {
  sample_name <- str_match(wang_matrices[i], "TNBC/.+(BC.+)_")
  sample_name <- sample_name[,2]
  
  count_matrix <- Matrix::readMM(wang_matrices[i])
  genes <- read.delim(wang_genes[[i]], header = FALSE)
  barcodes <- read.delim(wang_barcodes[[i]], header = FALSE)
  count_matrix <- as(as.matrix(count_matrix), "dgCMatrix")
  rownames(count_matrix) <- genes$V2
  colnames(count_matrix) <- barcodes$V1
  count_matrix <- count_matrix[!duplicated(genes$V2),]
  seurat_obj <- CreateSeuratObject(counts=count_matrix)
  
  seurat_obj$nUMI <- colSums(count_matrix)
  seurat_obj$nGene <- colSums(count_matrix != 0)
  seurat_obj$mitoRatio <- PercentageFeatureSet(object = seurat_obj, pattern = "^MT-")
  seurat_obj$mitoRatio <- seurat_obj@meta.data$mitoRatio / 100
  print(seurat_obj$nUMI %>% fivenum)
  print(seurat_obj$nGene %>% fivenum)
  print(seurat_obj$mitoRatio %>% fivenum)
  seurat_obj$log1p_counts <- log1p(seurat_obj$nCount_RNA)
  seurat_obj$log1p_n_gene <- log1p(seurat_obj$nFeature_RNA)
  seurat_obj$donor <- sample_name
  # Seems like nUMI was already filtered
  seurat_obj$low_quality <- mad_filter(seurat_obj, "log1p_counts", n_mads = 3) | mad_filter(seurat_obj, "log1p_n_gene", n_mads = 3) | mad_filter(seurat_obj, "mitoRatio", 3)
  # # According to their methods section
  # seurat_obj <- subset(x = seurat_obj,
  #                      subset= (nGene > 200) &
  #                        (nGene < 2500)
  #                        (mitoRatio < 0.2))
  seurat_obj <- subset(x = seurat_obj, subset = (low_quality != TRUE))
  seurat_obj <- Seurat::NormalizeData(seurat_obj)
  
  sce_obj <- Seurat::as.SingleCellExperiment(seurat_obj)
  sce_obj <- scDblFinder::scDblFinder(sce_obj)
  saveRDS(sce_obj, file=file.path(DATA_PATH, paste0("wang_2024/processed/TNBC/", sample_name, "_sce.rds")))
}
processed_paths <- Sys.glob(file.path(DATA_PATH, "wang_2024/processed/TNBC/*_sce.rds"))
wang_rds <- lapply(processed_paths, readRDS)
for (i in seq_along(wang_rds)) {
  rowData(wang_rds[[i]])$scDblFinder.selected <- NULL
}
wang_common_features <- purrr::reduce(lapply(wang_rds, rownames), intersect)
wang_rds <- lapply(wang_rds, function(x) x[rownames(x) %in% wang_common_features,])
wang_combined <- purrr::reduce(wang_rds, SingleCellExperiment::cbind)
wang_combined$batch <- "wang_2024"
wang_combined$celltype_major <- NA
wang_combined$celltype_minor <- NA
saveRDS(wang_combined, file=file.path(DATA_PATH, "wang_2024/processed/TNBC/combined.rds"))
