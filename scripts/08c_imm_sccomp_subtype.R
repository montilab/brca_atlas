library(tidyverse)
library(sccomp)

PATH <- file.path(Sys.getenv("MLAB"), "projects/brcameta/brca_atlas/")

# Loading Data
imm_seurat <- readRDS(file.path(PATH, "data/sc/imm_rpca_clean.rds"))
imm_metadata <- imm_seurat@meta.data
imm_metadata$sample_id <- rownames(imm_metadata)

# Adding Pam 50 data
pam50 <- read.csv(file.path(PATH, "data/metadata/pam50.csv"))[,c("donor", "pam50")]
imm_metadata$donor <- str_replace_all(imm_metadata$donor, "-", "_")
imm_metadata <- dplyr::left_join(x = imm_metadata, y = pam50, by = "donor")

# Binarizing factors
imm_metadata$sub_tnbc <- with(imm_metadata, dplyr::if_else(subtype_new == "TNBC", TRUE, FALSE))
imm_metadata$sub_er <- with(imm_metadata, dplyr::if_else(subtype_new == "ER+", TRUE, FALSE))
imm_metadata$sub_her2 <- with(imm_metadata, dplyr::if_else(subtype_new == "HER2+", TRUE, FALSE))
imm_metadata$pam_luma <- with(imm_metadata, dplyr::if_else(pam50 == "LumA", TRUE, FALSE))
imm_metadata$pam_lumb <- with(imm_metadata, dplyr::if_else(pam50 == "LumB", TRUE, FALSE))
imm_metadata$pam_norm <- with(imm_metadata, dplyr::if_else(pam50 == "Normal", TRUE, FALSE))
imm_metadata$pam_her2 <- with(imm_metadata, dplyr::if_else(pam50 == "Her2", TRUE, FALSE))
imm_metadata$pam_basal <- with(imm_metadata, dplyr::if_else(pam50 == "Basal", TRUE, FALSE))

all.equal(colnames(imm_seurat), imm_metadata$sample_id)
rownames(imm_metadata) <- imm_metadata$sample_id
imm_seurat@meta.data <- imm_metadata

# sccomp needs complete cases
luma_samples <- imm_metadata$sample_id[complete.cases(imm_metadata %>% 
                                                        dplyr::select(sample_id, pam_luma, cluster_annot, grade, batch, age))]
lumb_samples <- imm_metadata$sample_id[complete.cases(imm_metadata %>% 
                                                        dplyr::select(sample_id, pam_lumb, cluster_annot, grade, batch, age))]
norm_samples <- imm_metadata$sample_id[complete.cases(imm_metadata %>% 
                                                        dplyr::select(sample_id, pam_norm, cluster_annot, grade, batch, age))]
her2_samples <- imm_metadata$sample_id[complete.cases(imm_metadata %>% 
                                                        dplyr::select(sample_id, pam_her2, cluster_annot, grade, batch, age))]
basal_samples <- imm_metadata$sample_id[complete.cases(imm_metadata %>% 
                                                         dplyr::select(sample_id, pam_basal, cluster_annot, grade, batch, age))]
sub_her2_samples <- imm_metadata$sample_id[complete.cases(imm_metadata %>% 
                                                            dplyr::select(sample_id, sub_her2, cluster_annot, grade, batch, age))]
sub_er_samples <- imm_metadata$sample_id[complete.cases(imm_metadata %>% 
                                                          dplyr::select(sample_id, sub_er, cluster_annot, grade, batch, age))]
sub_tnbc_samples <- imm_metadata$sample_id[complete.cases(imm_metadata %>% 
                                                            dplyr::select(sample_id, sub_tnbc, cluster_annot, grade, batch, age))]

# 1. LumA 
imm_seurat_nona <- imm_seurat[, luma_samples]
print(imm_seurat_nona@meta.data$donor %>% unique %>% length)

sccomp_result =
  imm_seurat_nona |>
  sccomp_estimate(
    formula_composition = ~ pam_luma + grade + age + (1|batch),
    .sample =  donor,
    .cell_group = cluster_annot,
    bimodal_mean_variability_association = TRUE,
    cores = 16
  ) |>
  sccomp_remove_outliers(cores = 1) |> # Optional
  sccomp_test()

saveRDS(sccomp_result, file.path(PATH, "results/sccomp/imm/pam_luma.rds"))
print("Done with LumA Regression")

# 2. LumB
imm_seurat_nona <- imm_seurat[, lumb_samples]
print(imm_seurat_nona@meta.data$donor %>% unique %>% length)

sccomp_result =
  imm_seurat_nona |>
  sccomp_estimate(
    formula_composition = ~ pam_lumb + grade + age + (1|batch),
    .sample =  donor,
    .cell_group = cluster_annot,
    bimodal_mean_variability_association = TRUE,
    cores = 16
  ) |>
  sccomp_remove_outliers(cores = 1) |> # Optional
  sccomp_test()

saveRDS(sccomp_result, file.path(PATH, "results/sccomp/imm/pam_lumb.rds"))
print("Done with Lumb Regression")

# 3. Basal
imm_seurat_nona <- imm_seurat[, basal_samples]
print(imm_seurat_nona@meta.data$donor %>% unique %>% length)

sccomp_result =
  imm_seurat_nona |>
  sccomp_estimate(
    formula_composition = ~ pam_basal + grade + age + (1|batch),
    .sample =  donor,
    .cell_group = cluster_annot,
    bimodal_mean_variability_association = TRUE,
    cores = 16
  ) |>
  sccomp_remove_outliers(cores = 1) |> # Optional
  sccomp_test()

saveRDS(sccomp_result, file.path(PATH, "results/sccomp/imm/pam_basal.rds"))
print("Done with basal Regression")

# 4. Pam Her2 
imm_seurat_nona <- imm_seurat[, her2_samples]
print(imm_seurat_nona@meta.data$donor %>% unique %>% length)

sccomp_result =
  imm_seurat_nona |>
  sccomp_estimate(
    formula_composition = ~ pam_her2 + grade + age + (1|batch),
    .sample =  donor,
    .cell_group = cluster_annot,
    bimodal_mean_variability_association = TRUE,
    cores = 16
  ) |>
  sccomp_remove_outliers(cores = 1) |> # Optional
  sccomp_test()

saveRDS(sccomp_result, file.path(PATH, "results/sccomp/imm/pam_her2.rds"))
print("Done with Her2 Regression")

# 5. Norm
imm_seurat_nona <- imm_seurat[, norm_samples]
print(imm_seurat_nona@meta.data$donor %>% unique %>% length)

sccomp_result =
  imm_seurat_nona |>
  sccomp_estimate(
    formula_composition = ~ pam_norm + grade + age + (1|batch),
    .sample =  donor,
    .cell_group = cluster_annot,
    bimodal_mean_variability_association = TRUE,
    cores = 16
  ) |>
  sccomp_remove_outliers(cores = 1) |> # Optional
  sccomp_test()

saveRDS(sccomp_result, file.path(PATH, "results/sccomp/imm/pam_norm.rds"))
print("Done with Norm Regression")

# 6. TNBC
imm_seurat_nona <- imm_seurat[, sub_tnbc_samples]
print(imm_seurat_nona@meta.data$donor %>% unique %>% length)

sccomp_result =
  imm_seurat_nona |>
  sccomp_estimate(
    formula_composition = ~ sub_tnbc + grade + age + (1|batch),
    .sample =  donor,
    .cell_group = cluster_annot,
    bimodal_mean_variability_association = TRUE,
    cores = 16
  ) |>
  sccomp_remove_outliers(cores = 1) |> # Optional
  sccomp_test()

saveRDS(sccomp_result, file.path(PATH, "results/sccomp/imm/sub_tnbc.rds"))
print("Done with TNBC Regression")

# 7. Her2+
imm_seurat_nona <- imm_seurat[, sub_her2_samples]
print(imm_seurat_nona@meta.data$donor %>% unique %>% length)

sccomp_result =
  imm_seurat_nona |>
  sccomp_estimate(
    formula_composition = ~ sub_her2 + grade + age + (1|batch),
    .sample =  donor,
    .cell_group = cluster_annot,
    bimodal_mean_variability_association = TRUE,
    cores = 16
  ) |>
  sccomp_remove_outliers(cores = 1) |> # Optional
  sccomp_test()

saveRDS(sccomp_result, file.path(PATH, "results/sccomp/imm/sub_her2.rds"))
print("Done with Her2+ Regression")

# 8. ER+
imm_seurat_nona <- imm_seurat[, sub_er_samples]
print(imm_seurat_nona@meta.data$donor %>% unique %>% length)

sccomp_result =
  imm_seurat_nona |>
  sccomp_estimate(
    formula_composition = ~ sub_er + grade + age + (1|batch),
    .sample =  donor,
    .cell_group = cluster_annot,
    bimodal_mean_variability_association = TRUE,
    cores = 16
  ) |>
  sccomp_remove_outliers(cores = 1) |> # Optional
  sccomp_test()

saveRDS(sccomp_result, file.path(PATH, "results/sccomp/imm/sub_er.rds"))
print("Done with ER+ Regression")
