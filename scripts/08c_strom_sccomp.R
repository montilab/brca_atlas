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

# Adding metabolic data
all_metadata <- read.csv(file.path(PATH, "data/metadata/all_metadata.csv"), row.names = 1)
all_metadata <- all_metadata %>% dplyr::select(donor, obesity_prob, hardy_diff)
all_metadata$donor <- str_replace_all(all_metadata$donor, "-", "_")
strom_metadata <- dplyr::left_join(strom_metadata, all_metadata, by = "donor", multiple = "any")

# Binarizing factors
strom_metadata$grade_bin <- with(strom_metadata, dplyr::case_when(grade == 3 ~ "high",
                                                              grade == 2 ~ "low",
                                                              grade == 1 ~ "low"))
strom_metadata$ob_quant <- with(strom_metadata, cut(obesity_prob, quantile(obesity_prob, na.rm=TRUE),include.lowest=TRUE,labels=FALSE))
strom_metadata$ob_bin <- with(strom_metadata, if_else(ob_quant == 4, "Obese", "Non-Obese"))
strom_metadata$ir_quant <- with(strom_metadata, cut(hardy_diff, quantile(hardy_diff, na.rm=TRUE),include.lowest=TRUE,labels=FALSE))
strom_metadata$ir_bin <- with(strom_metadata, if_else(ir_quant == 4, "IR", "IS"))

strom_metadata$grade_bin <- factor(strom_metadata$grade_bin, levels = c("low", "high"))
strom_metadata$ob_bin <- factor(strom_metadata$ob_bin, levels = c("Non-Obese", "Obese"))
strom_metadata$ir_bin <- factor(strom_metadata$ir_bin, levels = c("IS", "IR"))

all.equal(colnames(strom_seurat), strom_metadata$sample_id)
rownames(strom_metadata) <- strom_metadata$sample_id
strom_seurat@meta.data <- strom_metadata

# sccomp needs complete cases
df_1 <- strom_metadata %>% dplyr::select(sample_id, cluster_annot, grade_bin, batch, subtype_new, age)
df_2 <- strom_metadata %>% dplyr::select(sample_id, cluster_annot, grade_bin, batch, pam50, age)
df_3 <- strom_metadata %>% dplyr::select(sample_id, cluster_annot, ob_bin, batch, subtype_new, age)
df_4 <- strom_metadata %>% dplyr::select(sample_id, cluster_annot, ob_bin, batch, pam50, age)
df_5 <- strom_metadata %>% dplyr::select(sample_id, cluster_annot, ir_bin, batch, subtype_new, age)
df_6 <- strom_metadata %>% dplyr::select(sample_id, cluster_annot, ir_bin, batch, pam50, age)

grade_subtype_df <- df_1[complete.cases(df_1),]
grade_pam_df <- df_2[complete.cases(df_2),]
ob_subtype_df <- df_3[complete.cases(df_3),]
ob_pam_df <- df_4[complete.cases(df_4),]
ir_subtype_df <- df_5[complete.cases(df_5),]
ir_pam_df <- df_6[complete.cases(df_6),]

grade_subtype_samples <- df_1$sample_id[complete.cases(df_1)]
grade_pam50_samples <- df_2$sample_id[complete.cases(df_2)]
ir_subtype_samples <- df_3$sample_id[complete.cases(df_3)]
ir_pam50_samples <- df_4$sample_id[complete.cases(df_4)]
ob_subtype_samples <- df_5$sample_id[complete.cases(df_5)]
ob_pam50_samples <- df_6$sample_id[complete.cases(df_6)]

# 1. Grade + Subtype Regression
strom_seurat_nona <- strom_seurat[, grade_subtype_samples]
print(strom_seurat_nona@meta.data$donor %>% unique %>% length)

sccomp_result = 
  strom_seurat_nona |>
  sccomp_estimate( 
    formula_composition = ~ grade_bin + subtype_new + age + (1|batch), 
    .sample =  donor, 
    .cell_group = cluster_annot, 
    bimodal_mean_variability_association = TRUE,
    cores = 1
  ) |> 
  sccomp_remove_outliers(cores = 1) |> # Optional
  sccomp_test()

saveRDS(sccomp_result, file.path(PATH, "results/sccomp/strom_grade_subtype.rds"))
print("Done with Grade/Subtype Regression")

# 2. Grade + pam50 Regression
strom_seurat_nona <- strom_seurat[, grade_pam50_samples]
print(strom_seurat_nona@meta.data$donor %>% unique %>% length)
sccomp_result = 
  strom_seurat_nona |>
  sccomp_estimate( 
    formula_composition = ~ grade_bin + pam50 + age + (1|batch), 
    .sample =  donor, 
    .cell_group = cluster_annot, 
    bimodal_mean_variability_association = TRUE,
    cores = 16 
  ) |> 
  sccomp_remove_outliers(cores = 1) |> # Optional
  sccomp_test()
saveRDS(sccomp_result, file.path(PATH, "results/sccomp/strom_grade_pam50.rds"))
print("Done with Grade/pam50 Regression")

# 3. IR/IS + Subtype Regression
strom_seurat_nona <- strom_seurat[, ir_subtype_samples]
print(strom_seurat_nona@meta.data$donor %>% unique %>% length)
sccomp_result = 
  strom_seurat_nona |>
  sccomp_estimate( 
    formula_composition = ~ ir_bin + subtype_new + age + (1|batch), 
    .sample =  donor, 
    .cell_group = cluster_annot, 
    bimodal_mean_variability_association = TRUE,
    cores = 16 
  ) |> 
  sccomp_remove_outliers(cores = 1) |> # Optional
  sccomp_test()
saveRDS(sccomp_result, file.path(PATH, "results/sccomp/strom_ir_subtype.rds"))
print("Done with IR/Subtype Regression")

# 4. IR/IS + pam Regression
strom_seurat_nona <- strom_seurat[, ir_pam50_samples]
print(strom_seurat_nona@meta.data$donor %>% unique %>% length)

sccomp_result = 
  strom_seurat_nona |>
  sccomp_estimate( 
    formula_composition = ~ ir_bin + pam50 + age + (1|batch), 
    .sample =  donor, 
    .cell_group = cluster_annot, 
    bimodal_mean_variability_association = TRUE,
    cores = 16 
  ) |> 
  sccomp_remove_outliers(cores = 1) |> # Optional
  sccomp_test()
saveRDS(sccomp_result, file.path(PATH, "results/sccomp/strom_ir_pam50.rds"))
print("Done with IR/pam50 Regression")

# 5. OB + subtype Regression
strom_seurat_nona <- strom_seurat[, ob_subtype_samples]
print(strom_seurat_nona@meta.data$donor %>% unique %>% length)

sccomp_result = 
  strom_seurat_nona |>
  sccomp_estimate( 
    formula_composition = ~ ob_bin + subtype_new + age + (1|batch), 
    .sample =  donor, 
    .cell_group = cluster_annot, 
    bimodal_mean_variability_association = TRUE,
    cores = 16 
  ) |> 
  sccomp_remove_outliers(cores = 1) |> # Optional
  sccomp_test()
saveRDS(sccomp_result, file.path(PATH, "results/sccomp/strom_ob_subtype.rds"))
print("Done with OB/Subtype Regression")

# 6. OB + pam50 Regression
strom_seurat_nona <- strom_seurat[, ob_pam50_samples]
print(strom_seurat_nona@meta.data$donor %>% unique %>% length)
sccomp_result = 
  strom_seurat_nona |>
  sccomp_estimate( 
    formula_composition = ~ ob_bin + pam50 + age + (1|batch), 
    .sample =  donor, 
    .cell_group = cluster_annot, 
    bimodal_mean_variability_association = TRUE,
    cores = 16 
  ) |> 
  sccomp_remove_outliers(cores = 1) |> # Optional
  sccomp_test()
saveRDS(sccomp_result, file.path(PATH, "results/sccomp/strom_ob_pam50.rds"))
print("Done with OB/pam50 Regression")