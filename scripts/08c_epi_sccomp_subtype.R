library(tidyverse)
library(sccomp)

PATH <- file.path(Sys.getenv("MLAB"), "projects/brcameta/brca_atlas/")

# Loading Data
epi_seurat <- readRDS(file.path(PATH, "data/sc/epi_rpca_subset.rds"))
epi_seurat <- subset(epi_seurat,
                     subset = RNA_snn_res.0.1 %in% setdiff(0:22, c(9,19,12,21,7,10,11)))
epi_metadata <- epi_seurat@meta.data
epi_metadata$sample_id <- rownames(epi_metadata)

# Adding Pam 50 data
pam50 <- read.csv(file.path(PATH, "data/metadata/pam50.csv"))[,c("donor", "pam50")]
epi_metadata$donor <- str_replace_all(epi_metadata$donor, "-", "_")
epi_metadata <- dplyr::left_join(x = epi_metadata, y = pam50, by = "donor")

# Binarizing factors
epi_metadata$sub_tnbc <- with(epi_metadata, dplyr::if_else(subtype_new == "TNBC", TRUE, FALSE))
epi_metadata$sub_er <- with(epi_metadata, dplyr::if_else(subtype_new == "ER+", TRUE, FALSE))
epi_metadata$sub_her2 <- with(epi_metadata, dplyr::if_else(subtype_new == "HER2+", TRUE, FALSE))
epi_metadata$pam_luma <- with(epi_metadata, dplyr::if_else(pam50 == "LumA", TRUE, FALSE))
epi_metadata$pam_lumb <- with(epi_metadata, dplyr::if_else(pam50 == "LumB", TRUE, FALSE))
epi_metadata$pam_norm <- with(epi_metadata, dplyr::if_else(pam50 == "Normal", TRUE, FALSE))
epi_metadata$pam_her2 <- with(epi_metadata, dplyr::if_else(pam50 == "Her2", TRUE, FALSE))
epi_metadata$pam_basal <- with(epi_metadata, dplyr::if_else(pam50 == "Basal", TRUE, FALSE))

all.equal(colnames(epi_seurat), epi_metadata$sample_id)
rownames(epi_metadata) <- epi_metadata$sample_id
epi_seurat@meta.data <- epi_metadata

# sccomp needs complete cases
luma_samples <- epi_metadata$sample_id[complete.cases(epi_metadata %>%
                                                        dplyr::select(sample_id, pam_luma, `RNA_snn_res.0.1`, grade, batch, age))]
lumb_samples <- epi_metadata$sample_id[complete.cases(epi_metadata %>%
                                                        dplyr::select(sample_id, pam_lumb, `RNA_snn_res.0.1`, grade, batch, age))]
norm_samples <- epi_metadata$sample_id[complete.cases(epi_metadata %>%
                                                        dplyr::select(sample_id, pam_norm, `RNA_snn_res.0.1`, grade, batch, age))]
her2_samples <- epi_metadata$sample_id[complete.cases(epi_metadata %>%
                                                        dplyr::select(sample_id, pam_her2, `RNA_snn_res.0.1`, grade, batch, age))]
basal_samples <- epi_metadata$sample_id[complete.cases(epi_metadata %>%
                                                         dplyr::select(sample_id, pam_basal, `RNA_snn_res.0.1`, grade, batch, age))]
sub_her2_samples <- epi_metadata$sample_id[complete.cases(epi_metadata %>%
                                                            dplyr::select(sample_id, sub_her2, `RNA_snn_res.0.1`, grade, batch, age))]
sub_er_samples <- epi_metadata$sample_id[complete.cases(epi_metadata %>%
                                                          dplyr::select(sample_id, sub_er, `RNA_snn_res.0.1`, grade, batch, age))]
sub_tnbc_samples <- epi_metadata$sample_id[complete.cases(epi_metadata %>%
                                                            dplyr::select(sample_id, sub_tnbc, `RNA_snn_res.0.1`, grade, batch, age))]

# 1. LumA
epi_seurat_nona <- epi_seurat[, luma_samples]
print(epi_seurat_nona@meta.data$donor %>% unique %>% length)

sccomp_result =
  epi_seurat_nona |>
  sccomp_estimate(
    formula_composition = ~ pam_luma + grade + age + (1|batch),
    .sample =  donor,
    .cell_group = `RNA_snn_res.0.1`,
    bimodal_mean_variability_association = TRUE,
    cores = 16
  ) |>
  sccomp_remove_outliers(cores = 1) |> # Optional
  sccomp_test()

saveRDS(sccomp_result, file.path(PATH, "results/sccomp/epi/pam_luma.rds"))
print("Done with LumA Regression")

# 2. LumB
epi_seurat_nona <- epi_seurat[, lumb_samples]
print(epi_seurat_nona@meta.data$donor %>% unique %>% length)

sccomp_result =
  epi_seurat_nona |>
  sccomp_estimate(
    formula_composition = ~ pam_lumb + grade + age + (1|batch),
    .sample =  donor,
    .cell_group = `RNA_snn_res.0.1`,
    bimodal_mean_variability_association = TRUE,
    cores = 16
  ) |>
  sccomp_remove_outliers(cores = 1) |> # Optional
  sccomp_test()

saveRDS(sccomp_result, file.path(PATH, "results/sccomp/epi/pam_lumb.rds"))
print("Done with Lumb Regression")

# 3. Basal
epi_seurat_nona <- epi_seurat[, basal_samples]
print(epi_seurat_nona@meta.data$donor %>% unique %>% length)

sccomp_result =
  epi_seurat_nona |>
  sccomp_estimate(
    formula_composition = ~ pam_basal + grade + age + (1|batch),
    .sample =  donor,
    .cell_group = `RNA_snn_res.0.1`,
    bimodal_mean_variability_association = TRUE,
    cores = 16
  ) |>
  sccomp_remove_outliers(cores = 1) |> # Optional
  sccomp_test()

saveRDS(sccomp_result, file.path(PATH, "results/sccomp/epi/pam_basal.rds"))
print("Done with basal Regression")

# 4. Pam Her2
epi_seurat_nona <- epi_seurat[, her2_samples]
print(epi_seurat_nona@meta.data$donor %>% unique %>% length)

sccomp_result =
  epi_seurat_nona |>
  sccomp_estimate(
    formula_composition = ~ pam_her2 + grade + age + (1|batch),
    .sample =  donor,
    .cell_group = `RNA_snn_res.0.1`,
    bimodal_mean_variability_association = TRUE,
    cores = 16
  ) |>
  sccomp_remove_outliers(cores = 1) |> # Optional
  sccomp_test()

saveRDS(sccomp_result, file.path(PATH, "results/sccomp/epi/pam_her2.rds"))
print("Done with Her2 Regression")

# 5. Norm
epi_seurat_nona <- epi_seurat[, norm_samples]
print(epi_seurat_nona@meta.data$donor %>% unique %>% length)

sccomp_result =
  epi_seurat_nona |>
  sccomp_estimate(
    formula_composition = ~ pam_norm + grade + age + (1|batch),
    .sample =  donor,
    .cell_group = `RNA_snn_res.0.1`,
    bimodal_mean_variability_association = TRUE,
    cores = 16
  ) |>
  sccomp_remove_outliers(cores = 1) |> # Optional
  sccomp_test()

saveRDS(sccomp_result, file.path(PATH, "results/sccomp/epi/pam_norm.rds"))
print("Done with Norm Regression")

# 6. TNBC
epi_seurat_nona <- epi_seurat[, sub_tnbc_samples]
print(epi_seurat_nona@meta.data$donor %>% unique %>% length)

sccomp_result =
  epi_seurat_nona |>
  sccomp_estimate(
    formula_composition = ~ sub_tnbc + grade + age + (1|batch),
    .sample =  donor,
    .cell_group = `RNA_snn_res.0.1`,
    bimodal_mean_variability_association = TRUE,
    cores = 16
  ) |>
  sccomp_remove_outliers(cores = 1) |> # Optional
  sccomp_test()

saveRDS(sccomp_result, file.path(PATH, "results/sccomp/epi/sub_tnbc.rds"))
print("Done with TNBC Regression")

# 7. Her2+
epi_seurat_nona <- epi_seurat[, sub_her2_samples]
print(epi_seurat_nona@meta.data$donor %>% unique %>% length)

sccomp_result =
  epi_seurat_nona |>
  sccomp_estimate(
    formula_composition = ~ sub_her2 + grade + age + (1|batch),
    .sample =  donor,
    .cell_group = `RNA_snn_res.0.1`,
    bimodal_mean_variability_association = TRUE,
    cores = 16
  ) |>
  sccomp_remove_outliers(cores = 1) |> # Optional
  sccomp_test()

saveRDS(sccomp_result, file.path(PATH, "results/sccomp/epi/sub_her2.rds"))
print("Done with Her2+ Regression")

# 8. ER+
epi_seurat_nona <- epi_seurat[, sub_er_samples]
print(epi_seurat_nona@meta.data$donor %>% unique %>% length)

sccomp_result =
  epi_seurat_nona |>
  sccomp_estimate(
    formula_composition = ~ sub_er + grade + age + (1|batch),
    .sample =  donor,
    .cell_group = `RNA_snn_res.0.1`,
    bimodal_mean_variability_association = TRUE,
    cores = 16
  ) |>
  sccomp_remove_outliers(cores = 1) |> # Optional
  sccomp_test()

saveRDS(sccomp_result, file.path(PATH, "results/sccomp/epi/sub_er.rds"))
print("Done with ER+ Regression")
