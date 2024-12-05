library(tidyverse)
library(Seurat)
library(sceasy)
library(SingleCellExperiment)
library(reticulate)
reticulate::use_condaenv("r-sceasy")
library(scMerge)
library(readxl)
options(Seurat.object.assay.version = "v5")
source("util.R")

DATA_PATH <- file.path(Sys.getenv("CBM"), "otherStudies/scRNAseq/brca")
PATH <- file.path(Sys.getenv("MLAB"), "projects/brcameta/brca_atlas")

# # 1. Reading SCE Data
wu_natgen_2021 <- readRDS(file.path(DATA_PATH, "wu_natgen_2021/wu_natgen_seurat.rds"))
wu_embo_2020 <- readRDS(file.path(DATA_PATH, "wu_embo_2020/wu_embo_seurat.rds"))
wu_genomemed_2021 <- readRDS(file.path(DATA_PATH, "wu_genomemed_2021/wu_genomemed_seurat.rds"))
pal_cancer <- readRDS(file.path(DATA_PATH, "pal_2021/pal_combined_seurat.rds"))
qian <- readRDS(file.path(DATA_PATH, "qian_2020/qian_seurat.rds"))
gao <- readRDS(file.path(DATA_PATH, "gao_2021/gao_seurat.rds"))
tietscher_combined <- readRDS(file.path(DATA_PATH, "tietscher_2023/tietscher_combined_seurat.rds"))
bassez_2021 <- readRDS(file.path(DATA_PATH, "bassez_2021/bassez_combined_seurat.rds"))
liu <- readRDS(file.path(DATA_PATH, "liu_2023/liu_seurat.rds"))
wang_combined <- readRDS(file.path(DATA_PATH, "wang_2024/wang_combined_seurat.rds"))

# Adding Batch Variable
wu_natgen_2021$batch <- "wu_natgen_2021"
wu_embo_2020$batch <- "wu_embo_2020"
wu_genomemed_2021$batch <- "wu_genomemed_2021"
pal_cancer$batch <- "pal_2021"
qian$batch <- "qian_2020"
gao$batch <- "gao_2021"
tietscher_combined$batch <- "tietscher_2023"
bassez_2021$batch <- "bassez_2021"
liu$batch <- "liu_2023"
wang_combined$batch <- "wang_2024"
wang_combined$subtype <- NA

# Adding subtype information
liu$subtype <- with(liu@meta.data, case_when(donor == "T2" ~ "HER2+",
                                                   donor == "T3" ~ "HER2+",
                                                   donor == "T6" ~ "HER2+",
                                                   donor == "T7" ~ "ER+"))

wang_combined$subtype <- with(wang_combined@meta.data, case_when(donor == "BC389_Tumor" ~ "ER+",
                                                     donor == "BC392_Tumor" ~ "ER+",
                                                     donor == "BC393_Tumor" ~ "ER+",
                                                     donor == "BC397_Tumor" ~ "ER+",
                                                     donor == "BC17086-12_Tumor" ~ "ER+",
                                                     donor == "BC17086-24_Tumor" ~ "TNBC",
                                                     donor == "BC17086-25_Tumor" ~ "TNBC",
                                                     donor == "BC17086-35_Tumor" ~ "TNBC",
                                                     donor == "BC258_Tumor" ~ "ER+",
                                                     donor == "BC302_Tumor" ~ "ER+",
                                                     donor == "BC394_Tumor" ~ "TNBC",
                                                     donor == "BC397_Tumor" ~ "ER+",
                                                     donor == "BC401_Tumor" ~ "ER+",
                                                     donor == "BC419_Tumor" ~ "TNBC",
                                                     donor == "BC428_Tumor" ~ "ER+"))

# Check if it's not working because of rowData
combined <- list(wu_natgen_2021, wu_embo_2020, wu_genomemed_2021,
                 pal_cancer, qian, gao, tietscher_combined, bassez_2021,
                 liu, wang_combined)

combined_rownames <- lapply(combined, rownames)

combined_intersect_rownames <- purrr::reduce(combined_rownames, intersect)

combined <- lapply(combined, function(x) x[rownames(x) %in% combined_intersect_rownames,])

for(i in seq_along(combined)) {
  combined[[i]]@meta.data <- combined[[i]]@meta.data[,c("celltype_major", "celltype_minor", "singler_pred", 
                                            "celltypist_pred", "subtype", "batch", "donor", 
                                            "scDblFinder.score", "scDblFinder.class")]
}

#combined_seurat <- purrr::reduce(combined, merge)
combined_seurat <- scCustomize::Merge_Seurat_List(combined, 
                                                  add.cell.ids = c("wu_natgen_", "wu_embo_", "wu_genomemed_",
                                                                   "pal_", "qian_", "gao_",
                                                                   "tietscher_", "bassez_",
                                                                   "liu_", "wang_"))


# May not need to remove this. Patient 4 has history of treatment.
combined_seurat <- combined_seurat[, combined_seurat$donor != "wu_embo_CID4513"]

# # Function to make values unique
# make_unique <- function(x) {
#   dups <- duplicated(x)
#   if (any(dups)) {
#     suffix <- seq_along(x)
#     x[dups] <- paste0(x[dups], "_", suffix[dups])
#   }
#   return(x)
# }
# 
# ids <- make_unique(colnames(combined_seurat))
# colnames(combined_seurat) <- ids
# rownames(combined_seurat@meta.data) <- ids
#saveRDS(combined_seurat, file.path(PATH, "data/sc/combined_seurat_filtered.rds"))

#combined_seurat <- readRDS(file.path(PATH, "data/sc/combined_seurat_filtered.rds"))
# Reconciling Cell Types
pdata_new <- combined_seurat@meta.data %>% as.data.frame
stopifnot(all.equal(rownames(pdata_new), colnames(combined_seurat)))

pdata_new$celltype_new <- author_new(pdata_new$celltype_major)

# Adding Broad Immune, Epithelial, Stromal labels
pdata_new$singleR_broad <- singler_broad(pdata_new$singler_pred)
pdata_new$author_broad <- author_broad(pdata_new$celltype_new)

# Adding more accurate subtype information
pdata_new$subtype <- with(pdata_new,
                          case_when(
                            (batch == "wu_genomemed_2021") & (donor == "wu_genomemed_CID4471") ~ "ER+",
                            (batch == "wu_genomemed_2021") & (donor == "wu_genomemed_CID44971") ~ "TNBC",
                            (batch == "wu_genomemed_2021") & (donor == "wu_genomemed_CID4513") ~ "TNBC",
                            .default = subtype
                          ))
pdata_new$subtype_new <- with(pdata_new,
                              case_when(str_detect(subtype, "ER\\+") ~ "ER+",
                                        str_detect(subtype, "HER2") ~ "HER2+",
                                        str_detect(subtype, "TNBC") ~ "TNBC",
                                        str_detect(subtype, "Triple") ~ "TNBC",
                                        subtype == "LumA" ~ "ER+",
                                        subtype == "LumB" ~ "ER+",
                                        subtype == "PR+ " ~ "ER+",
                                        .default = "Unassigned"))

# Adding more accurate age information
wang_metadata <- read_excel(file.path(PATH,"data/metadata/supplementary_data/wang_supplementary/mmc1.xlsx"))
wang_metadata$pid <- str_replace(wang_metadata$pid, "_", "-") 
wang_metadata$pid <- paste0(wang_metadata$pid, "_Tumor")
wang_metadata$pid[wang_metadata$pid == "BC1708625_Tumor"] <- "BC17086-25_Tumor"
wang_metadata$pid[wang_metadata$pid == "BC1708635_Tumor"] <- "BC17086-35_Tumor"
wang_metadata$pid[wang_metadata$pid == "BC1708624_Tumor"] <- "BC17086-24_Tumor"
wu_natgen_2021_age <- read_excel(file.path(PATH, "data/metadata/supplementary_data/wu_natgen_supplementary/41588_2021_911_MOESM4_ESM.xlsx"),
                                 sheet = "Supplementary Table 1", skip = 3)
wu_natgen_2021_age$donor_id <- paste0("wu_natgen_CID", wu_natgen_2021_age$`Case ID`)
wu_natgen_2021_age$donor_id <- str_replace(wu_natgen_2021_age$donor_id, pattern = "-1", replacement = "")
wu_natgen_2021_age$donor <- wu_natgen_2021_age$donor_id
wu_natgen_2021_age$age <- wu_natgen_2021_age$Age
wu_natgen_2021_age$grade <- wu_natgen_2021_age$Grade
wu_natgen_2021_age$stage <- wu_natgen_2021_age$Stage
wu_natgen_2021_age <- wu_natgen_2021_age[,c("donor", "age", "grade", "stage")]
pal_er <- read_excel(file.path(PATH, "data/metadata/supplementary_data/pal_supplementary/pal_age.xlsx"),
                     sheet = "ER")
pal_her2_tnbc <- read_excel(file.path(PATH, "data/metadata/supplementary_data/pal_supplementary/pal_age.xlsx"),
                            sheet = "HER2_TNBC")
pal_age <- data.frame(donor = c(pal_er$`Specimen ID`, pal_her2_tnbc$`Specimen ID`),
                      age = c(pal_er$`Patient Age`, pal_her2_tnbc$`Patient Age`),
                      grade = c(pal_er$Grade, pal_her2_tnbc$Grade),
                      stage = NA)
pal_age <- pal_age[!str_detect(pal_age$donor, "LN"),]
pal_age <- pal_age[!is.na(pal_age$donor),]
#pal_age <- drop_na(pal_age)
pal_age$donor <- str_extract(pal_age$donor, "(0\\d+.+)")
pal_age$donor <- str_remove(pal_age$donor,pattern = regex("-T"))
pal_age$donor <- paste0("pal_Patient ", pal_age$donor)
pal_age <- pal_age[,c("donor", "age", "grade", "stage")]
pal_wu_natgen <- rbind(pal_age, wu_natgen_2021_age)
tietscher <- read_excel(file.path(PATH, "data/metadata/supplementary_data/tietscher_supplementary/41467_2022_35238_MOESM4_ESM.xlsx"),
                        skip = 2)
tietscher$stage <- paste0(tietscher$`T Staging`, tietscher$`N Staging`, tietscher$`M Staging`)
tietscher$grade <- str_extract(tietscher$Grade, "\\d+") %>% as.numeric
tietscher$donor <- paste0("T", tietscher$`Patient ID`)
# range_to_midpoint <- function(data) {
#   # Split the range into two numbers
#   bounds <- as.numeric(unlist(strsplit(data, "-")))
#   # Calculate the midpoint
#   midpoint <- mean(bounds)
#   return(midpoint)
# }
# tietscher$age <- vapply(tietscher$`Age range at Surgery`, range_to_midpoint, numeric(1)) %>% unname
tietscher$age <- NA
tietscher <- tietscher %>% dplyr::select(c("donor", "age", "grade", "stage"))
wang_metadata$donor <- wang_metadata$pid
wang_metadata$age <- wang_metadata$Age
wang_metadata$stage <- wang_metadata$Stage
wang_metadata$grade <- wang_metadata$Grade
wang_metadata <- wang_metadata %>% dplyr::select(c("donor", "age", "grade", "stage"))
combined_supp_data <- rbind(pal_age, wu_natgen_2021_age, tietscher, wang_metadata)

pdata_new <- pdata_new %>% rownames_to_column(var="samplename")

pdata_new <- dplyr::left_join(x = pdata_new, y = combined_supp_data, by = "donor")
pdata_new$age <- with(pdata_new,
                              case_when(
                                (donor == "wu_genomemed_CID4471") ~ 55,
                                (donor == "wu_genomemed_CID44971") ~ 49,
                                (donor == "wu_genomemed_CID4513") ~ 73,
                                (donor == "wu_embo_CID44041") ~ 35, #P1
                                (donor == "wu_embo_CID44971") ~ 49, #P2
                                (donor == "wu_embo_CID44991") ~ 47, #P3
                                (donor == "wu_embo_CID4513") ~ 73, #P4
                                (donor == "wu_embo_CID4515") ~ 67, #P5
                                (donor == "GSM5457199") ~ 54,
                                (donor == "GSM5457202") ~ 56,
                                (donor == "GSM5457205") ~ 68,
                                (donor == "GSM5457208") ~ 44,
                                (donor == "GSM5457211") ~ 40,
                                (donor == "wu_natgen_CID4530N") ~ 42,
                                (donor == "wu_natgen_CID4290A") ~ 88,
                                (donor == "wu_natgen_CID44041") ~ 35,
                                (donor == "wu_natgen_CID44971") ~ 49,
                                (donor == "wu_natgen_CID44991") ~ 47,
                                (donor == "wu_natgen_CID45171") ~ 58,
                                (donor =="pal_Patient 0029-7C") ~ 59,
                                (donor =="pal_Patient 0029-9C") ~ 59,
                                (donor == "T2") ~ 36,
                                (donor == "T3") ~ 47,
                                (donor == "T6") ~ 46,
                                (donor == "T7") ~ 62,
                                .default = age
                              ))
pdata_new$grade <- pdata_new$grade %>% as.numeric
pdata_new$grade <- with(pdata_new,
                      case_when(
                        (donor == "wu_genomemed_CID4471") ~ 2,
                        (donor == "wu_genomemed_CID44971") ~ 3,
                        (donor == "wu_genomemed_CID4513") ~ 3,
                        (donor == "wu_embo_CID44041") ~ 3,
                        (donor == "wu_embo_CID44971") ~ 3,
                        (donor == "wu_embo_CID44991") ~ 3,
                        (donor == "wu_embo_CID4513") ~ 3,
                        (donor == "wu_embo_CID4515") ~ 3,
                        (donor == "T2") ~ 3,
                        (donor == "T3") ~ 3,
                        (donor == "T6") ~ 3,
                        (donor == "T7") ~ 3,
                        .default = grade
                      ))
pdata_new$stage <- with(pdata_new,
                        case_when(
                          (donor == "T2") ~ "IIB",
                          (donor == "T3") ~ "IIB",
                          (donor == "T6") ~ "IIIA",
                          (donor == "T7") ~ "IIB",
                          .default = stage
                        ))
pdata_new <- pdata_new %>% column_to_rownames(var="samplename")
age_data <- pdata_new %>% dplyr::select(age)
write.csv(age_data, file.path(PATH, "data/metadata/age.csv"), row.names = TRUE)
write.csv(pdata_new, file.path(PATH, "data/metadata/all_metadata.csv"), row.names = TRUE)


if (!all.equal(rownames(pdata_new), colnames(combined_seurat))) {
  rownames(pdata_new) <- colnames(combined_seurat)
}
combined_seurat@meta.data <- pdata_new

#saveRDS(combined_seurat, file=file.path(PATH, "data/combined_seurat.rds"))
saveRDS(combined_seurat_filtered, file.path(PATH, "data/sc/combined_seurat_filtered.rds"))
combined_seurat_filtered$celltypist_broad <- celltypist_broad(combined_seurat_filtered$celltypist_pred)
combined_seurat_filtered[["RNA"]] <- as(combined_seurat_filtered[["RNA"]], "Assay")

sceasy::convertFormat(combined_seurat_filtered, from="seurat", to="anndata", 
                      main_layer = "counts", transfer_layers="data",
                      outFile=file.path(PATH, "data/sc/combined_anndata_pp.h5ad"))
