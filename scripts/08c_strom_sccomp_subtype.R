library(tidyverse)
library(sccomp)

PATH <- file.path(Sys.getenv("MLAB"), "projects/brcameta/brca_atlas/")

# Loading Data
strom_seurat <- readRDS(file.path(PATH, "data/sc/strom_rpca_clean.rds"))
strom_metadata <- strom_seurat@meta.data
strom_metadata$sample_id <- rownames(strom_metadata)

# Adding Pam 50 data
pam50 <- read.csv(file.path(PATH, "data/metadata/pam50.csv"))[,c("donor", "pam50")]
strom_metadata$donor <- str_replace_all(strom_metadata$donor, "-", "_")
strom_metadata <- dplyr::left_join(x = strom_metadata, y = pam50, by = "donor")

# Binarizing factors
strom_metadata$sub_tnbc <- with(strom_metadata, dplyr::if_else(subtype_new == "TNBC", TRUE, FALSE))
strom_metadata$sub_er <- with(strom_metadata, dplyr::if_else(subtype_new == "ER+", TRUE, FALSE))
strom_metadata$sub_her2 <- with(strom_metadata, dplyr::if_else(subtype_new == "HER2+", TRUE, FALSE))
strom_metadata$pam_luma <- with(strom_metadata, dplyr::if_else(pam50 == "LumA", TRUE, FALSE))
strom_metadata$pam_lumb <- with(strom_metadata, dplyr::if_else(pam50 == "LumB", TRUE, FALSE))
strom_metadata$pam_norm <- with(strom_metadata, dplyr::if_else(pam50 == "Normal", TRUE, FALSE))
strom_metadata$pam_her2 <- with(strom_metadata, dplyr::if_else(pam50 == "Her2", TRUE, FALSE))
strom_metadata$pam_basal <- with(strom_metadata, dplyr::if_else(pam50 == "Basal", TRUE, FALSE))

all.equal(colnames(strom_seurat), strom_metadata$sample_id)
rownames(strom_metadata) <- strom_metadata$sample_id
strom_seurat@meta.data <- strom_metadata

# sccomp needs complete cases
luma_samples <- strom_metadata$sample_id[complete.cases(strom_metadata %>% 
                                                        dplyr::select(sample_id, pam_luma, cluster_annot, grade, batch, age))]
lumb_samples <- strom_metadata$sample_id[complete.cases(strom_metadata %>% 
                                                        dplyr::select(sample_id, pam_lumb, cluster_annot, grade, batch, age))]
norm_samples <- strom_metadata$sample_id[complete.cases(strom_metadata %>% 
                                                        dplyr::select(sample_id, pam_norm, cluster_annot, grade, batch, age))]
her2_samples <- strom_metadata$sample_id[complete.cases(strom_metadata %>% 
                                                        dplyr::select(sample_id, pam_her2, cluster_annot, grade, batch, age))]
basal_samples <- strom_metadata$sample_id[complete.cases(strom_metadata %>% 
                                                         dplyr::select(sample_id, pam_basal, cluster_annot, grade, batch, age))]
sub_her2_samples <- strom_metadata$sample_id[complete.cases(strom_metadata %>% 
                                                            dplyr::select(sample_id, sub_her2, cluster_annot, grade, batch, age))]
sub_er_samples <- strom_metadata$sample_id[complete.cases(strom_metadata %>% 
                                                          dplyr::select(sample_id, sub_er, cluster_annot, grade, batch, age))]
sub_tnbc_samples <- strom_metadata$sample_id[complete.cases(strom_metadata %>% 
                                                            dplyr::select(sample_id, sub_tnbc, cluster_annot, grade, batch, age))]

# 1. LumA 
strom_seurat_nona <- strom_seurat[, luma_samples]
print(strom_seurat_nona@meta.data$donor %>% unique %>% length)

sccomp_result =
  strom_seurat_nona |>
  sccomp_estimate(
    formula_composition = ~ pam_luma + grade + age + (1|batch),
    .sample =  donor,
    .cell_group = cluster_annot,
    bimodal_mean_variability_association = TRUE,
    cores = 16
  ) |>
  sccomp_remove_outliers(cores = 1) |> # Optional
  sccomp_test()

saveRDS(sccomp_result, file.path(PATH, "results/sccomp/strom/pam_luma.rds"))
print("Done with LumA Regression")

# 2. LumB
strom_seurat_nona <- strom_seurat[, lumb_samples]
print(strom_seurat_nona@meta.data$donor %>% unique %>% length)

sccomp_result =
  strom_seurat_nona |>
  sccomp_estimate(
    formula_composition = ~ pam_lumb + grade + age + (1|batch),
    .sample =  donor,
    .cell_group = cluster_annot,
    bimodal_mean_variability_association = TRUE,
    cores = 16
  ) |>
  sccomp_remove_outliers(cores = 1) |> # Optional
  sccomp_test()

saveRDS(sccomp_result, file.path(PATH, "results/sccomp/strom/pam_lumb.rds"))
print("Done with Lumb Regression")

# 3. Basal
strom_seurat_nona <- strom_seurat[, basal_samples]
print(strom_seurat_nona@meta.data$donor %>% unique %>% length)

sccomp_result =
  strom_seurat_nona |>
  sccomp_estimate(
    formula_composition = ~ pam_basal + grade + age + (1|batch),
    .sample =  donor,
    .cell_group = cluster_annot,
    bimodal_mean_variability_association = TRUE,
    cores = 16
  ) |>
  sccomp_remove_outliers(cores = 1) |> # Optional
  sccomp_test()

saveRDS(sccomp_result, file.path(PATH, "results/sccomp/strom/pam_basal.rds"))
print("Done with basal Regression")

# 4. Pam Her2 
strom_seurat_nona <- strom_seurat[, her2_samples]
print(strom_seurat_nona@meta.data$donor %>% unique %>% length)

sccomp_result =
  strom_seurat_nona |>
  sccomp_estimate(
    formula_composition = ~ pam_her2 + grade + age + (1|batch),
    .sample =  donor,
    .cell_group = cluster_annot,
    bimodal_mean_variability_association = TRUE,
    cores = 16
  ) |>
  sccomp_remove_outliers(cores = 1) |> # Optional
  sccomp_test()

saveRDS(sccomp_result, file.path(PATH, "results/sccomp/strom/pam_her2.rds"))
print("Done with Her2 Regression")

# 5. Norm
strom_seurat_nona <- strom_seurat[, norm_samples]
print(strom_seurat_nona@meta.data$donor %>% unique %>% length)

sccomp_result =
  strom_seurat_nona |>
  sccomp_estimate(
    formula_composition = ~ pam_norm + grade + age + (1|batch),
    .sample =  donor,
    .cell_group = cluster_annot,
    bimodal_mean_variability_association = TRUE,
    cores = 16
  ) |>
  sccomp_remove_outliers(cores = 1) |> # Optional
  sccomp_test()

saveRDS(sccomp_result, file.path(PATH, "results/sccomp/strom/pam_norm.rds"))
print("Done with Norm Regression")

# 6. TNBC
strom_seurat_nona <- strom_seurat[, sub_tnbc_samples]
print(strom_seurat_nona@meta.data$donor %>% unique %>% length)

sccomp_result =
  strom_seurat_nona |>
  sccomp_estimate(
    formula_composition = ~ sub_tnbc + grade + age + (1|batch),
    .sample =  donor,
    .cell_group = cluster_annot,
    bimodal_mean_variability_association = TRUE,
    cores = 16
  ) |>
  sccomp_remove_outliers(cores = 1) |> # Optional
  sccomp_test()

saveRDS(sccomp_result, file.path(PATH, "results/sccomp/strom/sub_tnbc.rds"))
print("Done with TNBC Regression")

# 7. Her2+
strom_seurat_nona <- strom_seurat[, sub_her2_samples]
print(strom_seurat_nona@meta.data$donor %>% unique %>% length)

sccomp_result =
  strom_seurat_nona |>
  sccomp_estimate(
    formula_composition = ~ sub_her2 + grade + age + (1|batch),
    .sample =  donor,
    .cell_group = cluster_annot,
    bimodal_mean_variability_association = TRUE,
    cores = 16
  ) |>
  sccomp_remove_outliers(cores = 1) |> # Optional
  sccomp_test()

saveRDS(sccomp_result, file.path(PATH, "results/sccomp/strom/sub_her2.rds"))
print("Done with Her2+ Regression")

# 8. ER+
strom_seurat_nona <- strom_seurat[, sub_er_samples]
print(strom_seurat_nona@meta.data$donor %>% unique %>% length)

sccomp_result =
  strom_seurat_nona |>
  sccomp_estimate(
    formula_composition = ~ sub_er + grade + age + (1|batch),
    .sample =  donor,
    .cell_group = cluster_annot,
    bimodal_mean_variability_association = TRUE,
    cores = 16
  ) |>
  sccomp_remove_outliers(cores = 1) |> # Optional
  sccomp_test()

saveRDS(sccomp_result, file.path(PATH, "results/sccomp/strom/sub_er.rds"))
print("Done with ER+ Regression")
