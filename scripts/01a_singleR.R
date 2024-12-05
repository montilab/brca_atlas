library(Seurat)
library(BiocParallel)
library(SingleR)
library(sceasy)
library(tidyverse)
library("AnnotationDbi")
library("org.Hs.eg.db")
library(parallel)
detectCores()
register(MulticoreParam(workers=16))

# library(doFuture)
# doFuture::registerDoFuture()
# future::plan("multisession", workers = 36)
# Setting globals to 100GB
# options(future.globals.maxSize = 100000 * 1024^2)

NAVIN_PATH <- file.path(Sys.getenv("CBM"), "otherStudies/scRNAseq/brca/navin_2023")

PATH <- file.path(Sys.getenv("MLAB"), "projects/brcameta/brca_atlas/")
cancer_paths <- read.delim(file.path(PATH, "scripts/01a_seurat_paths"), header = FALSE)$V1
id <- as.numeric(Sys.getenv("SGE_TASK_ID"))
path <- cancer_paths[[id]]

print(path)
dataset <- readRDS(path)

# Annotation with Navin 
navin_data <- readRDS(file.path(NAVIN_PATH, "local.rds"))
# Navin needs HGNC
hgncid <- mapIds(org.Hs.eg.db,
                  keys=rownames(navin_data),
                  column="SYMBOL",
                  keytype="ENSEMBL",
                  multiVals="first")
hgncid <- hgncid[!is.na(hgncid)]
ensembl_hgnc_df <- data.frame(hgnc = unname(hgncid), ensembl=names(hgncid))

rownames_filter <- rownames(navin_data) %in% ensembl_hgnc_df$ensembl
navin_data <- navin_data[rownames_filter,]
ref_mat <- navin_data@assays$RNA$data
rownames(ref_mat) <- ensembl_hgnc_df$hgnc

pred <- SingleR(test=dataset@assays@data$logcounts,
                ref=ref_mat,
                labels=navin_data$cell_type,
                de.method="wilcox",
                fine.tune=TRUE,
                BPPARAM=MulticoreParam(16))


dataset$singler_pred <- pred$pruned.labels

saveRDS(dataset, path)
print("Done with singleR annotation!")
#saveRDS(pred, file.path(PATH, "data/singleR_pred2.rds"))

