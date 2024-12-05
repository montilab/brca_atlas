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

# Adding metabolic data
all_metadata <- read.csv(file.path(PATH, "data/metadata/all_metadata.csv"), row.names = 1)
all_metadata <- all_metadata %>% dplyr::select(donor, obesity_prob, hardy_diff)
all_metadata$donor <- str_replace_all(all_metadata$donor, "-", "_")
imm_metadata <- dplyr::left_join(imm_metadata, all_metadata, by = "donor", multiple = "any")

# Binarizing factors
imm_metadata$grade_bin <- with(imm_metadata, dplyr::case_when(grade == 3 ~ "high",
                                                              grade == 2 ~ "low",
                                                              grade == 1 ~ "low"))
imm_metadata$ob_quant <- with(imm_metadata, cut(obesity_prob, quantile(obesity_prob, na.rm=TRUE),include.lowest=TRUE,labels=FALSE))
imm_metadata$ob_bin <- with(imm_metadata, if_else(ob_quant == 4, "Obese", "Non-Obese"))
imm_metadata$ir_quant <- with(imm_metadata, cut(hardy_diff, quantile(hardy_diff, na.rm=TRUE),include.lowest=TRUE,labels=FALSE))
imm_metadata$ir_bin <- with(imm_metadata, if_else(ir_quant == 4, "IR", "IS"))

imm_metadata$grade_bin <- factor(imm_metadata$grade_bin, levels = c("low", "high"))
imm_metadata$ob_bin <- factor(imm_metadata$ob_bin, levels = c("Non-Obese", "Obese"))
imm_metadata$ir_bin <- factor(imm_metadata$ir_bin, levels = c("IS", "IR"))

all.equal(colnames(imm_seurat), imm_metadata$sample_id)
rownames(imm_metadata) <- imm_metadata$sample_id
imm_seurat@meta.data <- imm_metadata

# sccomp needs complete cases
df_1 <- imm_metadata %>% dplyr::select(sample_id, cluster_annot, grade_bin, batch, subtype_new, age)
df_2 <- imm_metadata %>% dplyr::select(sample_id, cluster_annot, grade_bin, batch, pam50, age)
df_3 <- imm_metadata %>% dplyr::select(sample_id, cluster_annot, ob_bin, batch, subtype_new, age)
df_4 <- imm_metadata %>% dplyr::select(sample_id, cluster_annot, ob_bin, batch, pam50, age)
df_5 <- imm_metadata %>% dplyr::select(sample_id, cluster_annot, ir_bin, batch, subtype_new, age)
df_6 <- imm_metadata %>% dplyr::select(sample_id, cluster_annot, ir_bin, batch, pam50, age)

grade_subtype_samples <- df_1$sample_id[complete.cases(df_1)]
grade_pam50_samples <- df_2$sample_id[complete.cases(df_2)]
ir_subtype_samples <- df_3$sample_id[complete.cases(df_3)]
ir_pam50_samples <- df_4$sample_id[complete.cases(df_4)]
ob_subtype_samples <- df_5$sample_id[complete.cases(df_5)]
ob_pam50_samples <- df_6$sample_id[complete.cases(df_6)]

# 1. Grade + Subtype Regression
imm_seurat_nona <- imm_seurat[, grade_subtype_samples]
print(imm_seurat_nona@meta.data$donor %>% unique %>% length)

# sccomp_result =
#   imm_seurat_nona |>
#   sccomp_estimate(
#     formula_composition = ~ grade_bin + subtype_new + age + (1|batch),
#     .sample =  donor,
#     .cell_group = cluster_annot,
#     bimodal_mean_variability_association = TRUE,
#     cores = 16
#   ) |>
#   sccomp_remove_outliers(cores = 1) |> # Optional
#   sccomp_test()
# 
# saveRDS(sccomp_result, file.path(PATH, "results/sccomp/imm/grade_subtype.rds"))
# print("Done with Grade/Subtype Regression")

# 1. b Subsetting by batch
# for (dataset in unique(imm_seurat_nona$batch)) {
#   imm_seurat_batch <- subset(imm_seurat_nona, subset = (batch == dataset))
#   if (n_distinct(imm_seurat_batch$grade_bin) < 2 | n_distinct(imm_seurat_batch$pam50) < 2 | n_distinct(imm_seurat_batch$donor) <=3 ) {
#     print(paste0("Skipping ", dataset))
#     next
#   }
#   print(dataset)
#   print(imm_seurat_batch)
#   print(imm_seurat_batch@meta.data$donor %>% unique %>% length)
#   
#   sccomp_result =
#     imm_seurat_batch |>
#     sccomp_estimate(
#       formula_composition = ~ grade_bin + pam50 + age,
#       .sample =  donor,
#       .cell_group = cluster_annot,
#       bimodal_mean_variability_association = TRUE,
#       cores = 16
#     ) |>
#     sccomp_remove_outliers(cores = 1) |> # Optional
#     sccomp_test()
# 
#   saveRDS(sccomp_result, file.path(PATH, paste0("results/sccomp/imm/",dataset,"_pam50_subtype.rds")))
# }

# 2. Grade + pam50 Regression
imm_seurat_nona <- imm_seurat[, grade_pam50_samples]
print(imm_seurat_nona@meta.data$donor %>% unique %>% length)
sccomp_result =
  imm_seurat_nona |>
  sccomp_estimate(
    formula_composition = ~ grade_bin + pam50 + age + (1|batch),
    formula_variability = ~ grade_bin + pam50 + age + (1|batch),
    .sample =  donor,
    .cell_group = cluster_annot,
    bimodal_mean_variability_association = TRUE,
    cores = 16
  ) |>
  sccomp_remove_outliers(cores = 1) |> # Optional
  sccomp_test()
saveRDS(sccomp_result, file.path(PATH, "results/sccomp/imm/grade_pam50_v.rds"))
# saveRDS(sccomp_result, file.path(PATH, "results/sccomp/imm/grade_pam50.rds"))
print("Done with Grade/pam50 Regression")
# 
# # 3. IR/IS + Subtype Regression
# imm_seurat_nona <- imm_seurat[, ir_subtype_samples]
# print(imm_seurat_nona@meta.data$donor %>% unique %>% length)
# sccomp_result =
#   imm_seurat_nona |>
#   sccomp_estimate(
#     formula_composition = ~ ir_bin + subtype_new + age + (1|batch),
#     .sample =  donor,
#     .cell_group = cluster_annot,
#     bimodal_mean_variability_association = TRUE,
#     cores = 16
#   ) |>
#   sccomp_remove_outliers(cores = 1) |> # Optional
#   sccomp_test()
# saveRDS(sccomp_result, file.path(PATH, "results/sccomp/imm/ir_subtype.rds"))
# print("Done with IR/Subtype Regression")
# 
# # 4. IR/IS + pam Regression
# imm_seurat_nona <- imm_seurat[, ir_pam50_samples]
# print(imm_seurat_nona@meta.data$donor %>% unique %>% length)
# 
# sccomp_result =
#   imm_seurat_nona |>
#   sccomp_estimate(
#     formula_composition = ~ ir_bin + pam50 + age + (1|batch),
#     .sample =  donor,
#     .cell_group = cluster_annot,
#     bimodal_mean_variability_association = TRUE,
#     cores = 16
#   ) |>
#   sccomp_remove_outliers(cores = 1) |> # Optional
#   sccomp_test()
# saveRDS(sccomp_result, file.path(PATH, "results/sccomp/imm/ir_pam.rds"))
# print("Done with IR/pam50 Regression")
# 
# # 5. OB + subtype Regression
# imm_seurat_nona <- imm_seurat[, ob_subtype_samples]
# print(imm_seurat_nona@meta.data$donor %>% unique %>% length)
# 
# sccomp_result =
#   imm_seurat_nona |>
#   sccomp_estimate(
#     formula_composition = ~ ob_bin + subtype_new + age + (1|batch),
#     .sample =  donor,
#     .cell_group = cluster_annot,
#     bimodal_mean_variability_association = TRUE,
#     cores = 16
#   ) |>
#   sccomp_remove_outliers(cores = 1) |> # Optional
#   sccomp_test()
# saveRDS(sccomp_result, file.path(PATH, "results/sccomp/imm/ob_subtype.rds"))
# print("Done with OB/Subtype Regression")
# 
# # 6. OB + pam50 Regression
# imm_seurat_nona <- imm_seurat[, ob_pam50_samples]
# print(imm_seurat_nona@meta.data$donor %>% unique %>% length)
# sccomp_result =
#   imm_seurat_nona |>
#   sccomp_estimate(
#     formula_composition = ~ ob_bin + pam50 + age + (1|batch),
#     .sample =  donor,
#     .cell_group = cluster_annot,
#     bimodal_mean_variability_association = TRUE,
#     cores = 16
#   ) |>
#   sccomp_remove_outliers(cores = 1) |> # Optional
#   sccomp_test()
# saveRDS(sccomp_result, file.path(PATH, "results/sccomp/imm/ob_pam50.rds"))
# print("Done with OB/pam50 Regression")
