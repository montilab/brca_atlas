library(tidyverse)
library(Seurat)
library(scCustomize)
library(ggrepel)
library(AUCell)
library(hypeR)
library(CytoTRACE2)
library(openxlsx)
library(ComplexHeatmap)
library(circlize)


PATH <- file.path(Sys.getenv("MLAB"), "projects/brcameta/brca_atlas/")

epi_rpca <- readRDS(file.path(PATH, "data/sc/epi_rpca_subset.rds")) 

# 1. Scoring with Hallmark/Reactome
# hallmark_genesets <- hypeR::msigdb_download(species="Homo sapiens", category="H")
# reactome_genesets <- hypeR::enrichr_download("Reactome_2022" )
# 
#K2_sce <- as.SingleCellExperiment(epi_rpca)
# cells_auc_reactome <- AUCell_run(K2_sce, reactome_genesets)
# saveRDS(cells_auc_reactome, file.path(PATH, "data/auc/epi_reactome.rds"))
# 
# cells_auc_hallmark <- AUCell_run(K2_sce, hallmark_genesets)
# saveRDS(cells_auc_hallmark, file.path(PATH, "data/auc/epi_hallmark.rds"))
#
# 2. Scoring with Proliferation/Inflammation
# prolif <- read.delim(file.path(PATH, "data/signatures/K2_TILS/PMID25848820_inflammationBreastCancer.txt"), header = FALSE)
# inflam <- read.delim(file.path(PATH, "data/signatures/K2_TILS/PMID22028643_proliferationBreastCancer.txt"), header = FALSE)
# cells_auc_prolif_inflam <- AUCell_run(K2_sce, list(proliferation = prolif$V1, inflammation = inflam$V1))
# saveRDS(cells_auc_prolif_inflam, file.path(PATH, "data/auc/epi_prolif_inflam.rds"))
#
# 3. Scoring with EMT
# emt_tan_sigs <- readRDS(file.path(PATH, "data/signatures/cancer_epithelial/tan_2014/tan_sigs.rds"))
# emt_winkler_sigs <- readRDS(file.path(PATH, "data/signatures/cancer_epithelial/winkler_2024/emp_sigs.rds"))
# cells_auc_emt <- AUCell_run(K2_sce, c(emt_tan_sigs, emt_winkler_sigs))
# saveRDS(cells_auc_emt, file.path(PATH, "data/auc/epi_emt.rds"))
#
# 4. Scoring with Healthy Breast Markers
# epi_all <- readRDS(file.path(PATH, "brca_atlas_validation/data/sigs/all_markers.rds"))
# cells_auc_hbca <- AUCell_run(K2_sce, epi_all)
# saveRDS(cells_auc_hbca, file.path(PATH, "data/auc/epi_hbca.rds"))
#
# 5. Scoring with cytotrace
# expression_matrix <- GetAssayData(epi_rpca, slot = "counts")
# 
# cytotrace2_result <- cytotrace2(expression_matrix, 
#                                 species = "human",
#                                 ncores = 16)
# 
# saveRDS(cytotrace2_result, file.path(PATH, "data/cytotrace/cytotrace_results.rds"))

# 6. Adding scores to a spreadsheet
combined_malig_data <- readRDS(file.path(PATH, "data/infercnv/combined_infercnv_data.rds"))
combined_malig_data <- combined_malig_data[,c("cell_id", "score", "annotation", "outside_norm_90")]
cytotrace_data <- readRDS(file.path(PATH, "data/cytotrace/cytotrace_results.rds"))
cytotrace_data <- rownames_to_column(cytotrace_data, "cell_id")
emp_scores <- readRDS(file.path(PATH, "data/auc/epi_emt.rds"))
hbca_scores <- readRDS(file.path(PATH, "data/auc/epi_hbca.rds"))
hallmark_scores <- readRDS(file.path(PATH, "data/auc/epi_hallmark.rds"))
pam50 <- read.csv(file.path(PATH, "data/metadata/pam50.csv"))[,c("donor", "pam50")]

emp_df <- data.frame(epi_breast = emp_scores@assays@data$AUC["epi_up",] - emp_scores@assays@data$AUC["epi_dn",],
                     int_breast = emp_scores@assays@data$AUC["int_up",] - emp_scores@assays@data$AUC["int_dn",],
                     mes_breast = emp_scores@assays@data$AUC["mes_up",] - emp_scores@assays@data$AUC["mes_dn",])
emp_df$cell_id <- colnames(emp_scores)
hbca_df <- as.data.frame(t(hbca_scores@assays@data$AUC[1:11,]))
hbca_df <- rownames_to_column(hbca_df, "cell_id")
hallmark_df <- as.data.frame(t(hallmark_scores@assays@data$AUC))
hallmark_df <- rownames_to_column(hallmark_df, "cell_id")

epi_metadata <- epi_rpca@meta.data
epi_metadata <- rownames_to_column(epi_metadata, "cell_id")
epi_metadata$donor <- str_replace_all(epi_metadata$donor, "-", "_")
epi_metadata <- dplyr::left_join(x = epi_metadata, y = pam50, by = "donor")

all_df <- list(epi_metadata, combined_malig_data, cytotrace_data, emp_df, hbca_df, hallmark_df)
merged_df <- Reduce(function(x, y) dplyr::left_join(x, y, by = "cell_id"), all_df)
merged_df_filtered <- merged_df[, c(1,14,15,16,23,30,31,33:35,39:102)]
colnames(merged_df_filtered) <- clean_hallmark_names(colnames(merged_df_filtered))
colnames(merged_df_filtered)[c(2,5,7,8,9,10,11,12,13)] <- c("Subtype","Cluster","InferCNV_Score", 
                                                            "InferCNV_Malig", "Cytotrace_Score", "Cytotrace_Class",
                                                            "Epi. State", "Int. State", "Mes. State")
saveRDS(merged_df_filtered, file = file.path(PATH, "results/annotation/epi/scores_cell.rds"))

# Fisher Test
grade_freq <- epi_metadata$grade %>% table %>% prop.table
subtype_freq <- epi_metadata$subtype_new %>% table %>% prop.table
# 
# # Function to perform Fisher's test on a single row
# row_fisher_test <- function(row, alternative = "greater") {
#   # Create 2x2 matrix from the row
#   mat <- matrix(as.numeric(row), nrow = 2, byrow = TRUE)
#   
#   # Perform Fisher's exact test
#   result <- fisher.test(mat, alternative = alternative)
#   
#   # Return p-value (you can modify this to return other statistics if needed)
#   return(list(odds_ratio = result$estimate %>% unname, 
#               p_value = result$p.value))
# }
# 
# grade_count_df <- merged_df %>% dplyr::group_by(RNA_snn_res.0.2, grade) %>% dplyr::summarize(Count = n())
# grade_total_df <- merged_df %>% dplyr::group_by(RNA_snn_res.0.2) %>% dplyr::summarize(Total = n())
# grade_merged_df <- dplyr::left_join(grade_count_df, grade_total_df, by = c("RNA_snn_res.0.2"))
# grade_merged_df <- grade_merged_df %>% drop_na()
# grade_merged_df <- grade_merged_df %>% dplyr::mutate(pop_count = case_when(grade == 1 ~ grade_freq[[1]],
#                                                                            grade == 2 ~ grade_freq[[2]],
#                                                                            grade == 3 ~ grade_freq[[3]]),
#                                                      pop_total = sum(grade_freq))
# 
# fisher_grade_values <- apply(grade_merged_df[,3:6], 1, row_fisher_test)
# fisher_grade_values <- as.data.frame(rbindlist(fisher_grade_values))
# grade_merged_df <- cbind(grade_merged_df, fisher_grade_values)
# grade_merged_df <- grade_merged_df %>% 
#   dplyr::select("RNA_snn_res.0.2", "grade", "odds_ratio") %>% 
#   tidyr::pivot_wider(names_from = "grade",
#                      values_from = "odds_ratio")
# 
# 
# subtype_count_df <- merged_df %>% dplyr::group_by(RNA_snn_res.0.2, subtype_new) %>% dplyr::summarize(Count = n())
# subtype_total_df <- merged_df %>% dplyr::group_by(RNA_snn_res.0.2) %>% dplyr::summarize(Total = n())
# subtype_merged_df <- dplyr::left_join(subtype_count_df, subtype_total_df, by = c("RNA_snn_res.0.2"))
# subtype_merged_df <- subtype_merged_df %>% dplyr::filter(subtype_new != "Unassigned")
# subtype_merged_df <- subtype_merged_df %>% dplyr::mutate(pop_count = case_when(subtype_new == "ER+" ~ subtype_freq[["ER+"]],
#                                                                                subtype_new == "HER2+" ~ subtype_freq[["HER2+"]],
#                                                                                subtype_new == "TNBC" ~ subtype_freq[["TNBC"]]),
#                                                          pop_total = sum(subtype_freq))
# 
# fisher_subtype_values <- apply(subtype_merged_df[,3:6], 1, row_fisher_test)
# fisher_subtype_values <- as.data.frame(rbindlist(fisher_subtype_values))
# subtype_merged_df <- cbind(subtype_merged_df, fisher_subtype_values)
# subtype_merged_df <- subtype_merged_df %>% 
#   dplyr::select("RNA_snn_res.0.2", "subtype_new", "odds_ratio") %>% 
#   tidyr::pivot_wider(names_from = "subtype_new",
#                      values_from = "odds_ratio")

metadata_df <- merged_df %>%
  dplyr::group_by(RNA_snn_res.0.2) %>%
  dplyr::summarize(prop_grade1 = mean(grade == 1, na.rm = TRUE),
                   prop_grade2 = mean(grade == 2, na.rm = TRUE),
                   prop_grade3 = mean(grade == 3, na.rm = TRUE),
                   prop_er = mean(subtype_new == "ER+", na.rm = TRUE),
                   prop_her2 = mean(subtype_new == "HER2+", na.rm = TRUE),
                   prop_tnbc = mean(subtype_new == "TNBC", na.rm = TRUE),
                   mean_age = mean(age, na.rm = TRUE))

# grouped_dfs <- list(metadata_df, subtype_merged_df, grade_merged_df)
# grouped_dfs <- list(merged_df, metadata_df)
# metadata_df <- Reduce(function(x, y) dplyr::left_join(x, y, by = "RNA_snn_res.0.2"), grouped_dfs)


cyto_malig_df <- merged_df %>% dplyr::select(RNA_snn_res.0.2, 
                                             score, outside_norm_90, 
                                             CytoTRACE2_Score, CytoTRACE2_Potency, 
                                             epi_breast, int_breast, mes_breast)
cyto_malig_df <- cyto_malig_df %>% 
  dplyr::group_by(RNA_snn_res.0.2) %>% 
  dplyr::summarise(prop_malignant = mean(outside_norm_90, na.rm = TRUE),
                   mean_cytotrace = mean(CytoTRACE2_Score),
                   prop_diff = mean(CytoTRACE2_Potency == "Differentiated"),
                   prop_uni = mean(CytoTRACE2_Potency == "Unipotent"),
                   prop_oligo = mean(CytoTRACE2_Potency == "Oligopotent"),
                   prop_pluri = mean(CytoTRACE2_Potency == "Pluripotent"),
                   prop_toti = mean(CytoTRACE2_Potency == "Totipotent"),
                   mean_epi_state = mean(epi_breast),
                   mean_int_state = mean(int_breast),
                   mean_mes_state = mean(mes_breast))
hbca_hallmark_df <- merged_df[,c(23, 41:101)]
hbca_hallmark_df <- hbca_hallmark_df %>% 
  dplyr::group_by(RNA_snn_res.0.2) %>%
  dplyr::summarise_all(mean)

grouped_dfs <- list(metadata_df, cyto_malig_df, hbca_hallmark_df)
combined_grouped_df <- Reduce(function(x, y) dplyr::left_join(x, y, by = "RNA_snn_res.0.2"), grouped_dfs)
combined_grouped_df$RNA_snn_res.0.2 <- as.numeric(as.character(combined_grouped_df$RNA_snn_res.0.2))
combined_grouped_df <- combined_grouped_df %>% dplyr::arrange(RNA_snn_res.0.2)

# Create a new workbook
wb <- createWorkbook()

addWorksheet(wb, sheetName = "Cluster Scores")
writeData(wb, sheet = "Cluster Scores", x = combined_grouped_df)

# Save the workbook
saveWorkbook(wb, file.path(PATH, "results/annotation/epi_scores.xlsx"), overwrite = TRUE)
saveRDS(combined_grouped_df, file.path(PATH, "results/annotation/epi_scores.rds"))

# 1. Epi Scores Heatmap
# Complex Heatmaps of scores
heatmap_dat <- combined_grouped_df %>% 
  dplyr::select(-c(prop_pluri, prop_toti)) %>%
  tibble::column_to_rownames(var = "RNA_snn_res.0.2")
heatmap_dat[,-c(1:8, 10:12)] <- scale(heatmap_dat[,-c(1:8, 10:12)])
heatmap_dat <- heatmap_dat %>% as.matrix %>% t

heatmap_labels <- rownames(heatmap_dat)
heatmap_labels <- clean_hallmark_names(heatmap_labels)
heatmap_labels[1:15] <- c("Grade 1 (6%)", 
                         "Grade 2 (25%)",
                         "Grade 3 (67%)",
                         "ER+ (56%)",
                         "HER2+ (10%)",
                         "TNBC (28%)",
                         "Mean Age",
                         "Prop. Malignant",
                        "Mean Cytotrace",
                        "Prop. Diff",
                        "Prop. Uni",
                        "Prop. Oligo",
                        "Mean Epi. State",
                        "Mean Int. State", 
                        "Mean Mes. State")
rownames(heatmap_dat) <- heatmap_labels

# Define two color schemes
col_scores <- colorRamp2(c(-3, 0, 3), c("blue", "white", "red"))
col_props <- colorRamp2(c(0, 0.5, 1),  c("#E5F5E0", "#A1D99B", "#31A354"))
col_age <- colorRamp2(c(min(heatmap_dat["Mean Age",]), 
                        median(heatmap_dat["Mean Age",]), 
                        max(heatmap_dat["Mean Age",])), c(c("#FEE6CE", "#FDAE6B", "#E6550D")))

ht <- Heatmap(heatmap_dat,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        cell_fun = function(j, i, x, y, width, height, fill) {
          if (i %in% c(1:6, 8, 10:12)) {
            grid.rect(x, y, width, height, gp = gpar(fill = col_props(heatmap_dat[i,j]), col = NA))
          } else if (i == 7) {
            grid.rect(x, y, width, height, gp = gpar(fill = col_age(heatmap_dat[i, j]), col = NA))
          } else {
            grid.rect(x, y, width, height, gp = gpar(fill = col_scores(heatmap_dat[i, j]), col = NA))
          }
          grid.text(sprintf("%.1f", heatmap_dat[i, j]), x, y, gp = gpar(fontsize = 3))
          },
        row_split = c(rep("Patient", 7), 
                      "InferCNV", 
                      rep("Cytotrace", 4), 
                      rep("EMP", 3),
                      rep("HBCA", 11),
                      rep("Hallmark", 50)),
        row_gap = unit(2, "mm"),
        row_title_gp = gpar(fontsize = 5),
        row_names_gp = gpar(fontsize = 5),
        column_names_gp = gpar(fontsize = 5), 
        column_names_rot = 45,
        heatmap_legend_param = list(
          labels_gp = gpar(fontsize = 5),  # Adjust size of legend labels
          title_gp = gpar(fontsize = 5)    # Adjust size of legend title
        ),
        name = "Scores",
        border = TRUE,
        show_heatmap_legend = FALSE)
# Create custom legends
lgd_score <- Legend(col_fun = col_scores, title = "Scores")
lgd_prop <- Legend(col_fun = col_props, title = "Proportions")
lgd_age <- Legend(col_fun = col_age, title = "Age")

# Draw the heatmap with custom legends

png(file.path(PATH, "results/annotation/epi_scores_heatmap.png"), width=1500, height=2000, res = 300)
draw(ht, heatmap_legend_list = list(lgd_score, lgd_prop, lgd_age))
dev.off()

# 2. Subset + Clustered Heatmap
minmax_scale <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}
# Heatmap without last four clusters
heatmap_dat <- combined_grouped_df %>% 
  dplyr::select(-c(prop_pluri, prop_toti)) %>%
  tibble::column_to_rownames(var = "RNA_snn_res.0.2")
heatmap_dat <- heatmap_dat[-c(16:19),]
heatmap_dat[,-c(1:8, 10:12)] <- apply(heatmap_dat[,-c(1:8, 10:12)], 2, minmax_scale)
#heatmap_dat[,-c(1:8, 10:12)] <- scale(heatmap_dat[,-c(1:8, 10:12)])

heatmap_dat <- heatmap_dat %>% as.matrix %>% t

heatmap_labels <- rownames(heatmap_dat)
heatmap_labels <- str_remove(heatmap_labels, "HALLMARK_")
heatmap_labels <- str_replace_all(heatmap_labels, pattern = "_", replacement = " ")
heatmap_labels[27:76] <- str_to_sentence(heatmap_labels[27:76])
heatmap_labels[27:76] <- str_replace(heatmap_labels[27:76], "Pi3k" , "PI3K") %>%
  str_replace(., "Wnt", "WNT") %>%
  str_replace(., "Uv", "UV") %>%
  str_replace(., "Tgf beta", "TGF Beta") %>% 
  str_replace(., "Tnfa", "TNFA") %>%
  str_replace(., "akt", "Akt") %>%
  str_replace(., "Mtorc1", "MTORC1") %>%
  str_replace(., "stat", "STAT") %>%
  str_replace(., "Myc", "MYC") %>%
  str_replace(., "jak", "JAK") %>%
  str_replace(., "G2m", "G2M") %>%
  str_replace(., "Il", "IL") %>%
  str_replace(., "E2f", "E2F") %>%
  str_replace(., "Dna", "DNA") %>%
  str_replace(., "Kras", "KRAS")
heatmap_labels[1:15] <- c("Grade 1 (6%)", 
                          "Grade 2 (25%)",
                          "Grade 3 (67%)",
                          "ER+ (56%)",
                          "HER2+ (10%)",
                          "TNBC (28%)",
                          "Mean Age",
                          "Prop. Malignant",
                          "Mean Cytotrace",
                          "Prop. Diff",
                          "Prop. Uni",
                          "Prop. Oligo",
                          "Mean Epi. State",
                          "Mean Int. State", 
                          "Mean Mes. State")
rownames(heatmap_dat) <- heatmap_labels
saveRDS(heatmap_dat, file = file.path(PATH, "results/annotation/epi/heatmap_dat_subset.rds"))

# Define two color schemes
# col_scores <- colorRamp2(c(-3, 0, 3), c("blue", "white", "red"))
col_scores <- colorRamp2(c(0, 1), c("white", "red"))
col_props <- colorRamp2(c(0, 0.5, 1),  c("#E5F5E0", "#A1D99B", "#31A354"))
col_age <- colorRamp2(c(min(heatmap_dat["Mean Age",]), 
                        median(heatmap_dat["Mean Age",]), 
                        max(heatmap_dat["Mean Age",])), c(c("#FEE6CE", "#FDAE6B", "#E6550D")))

ht <- Heatmap(heatmap_dat,
              cluster_rows = FALSE,
              cluster_columns = TRUE,
              cell_fun = function(j, i, x, y, width, height, fill) {
                if (i %in% c(1:6, 8, 10:12)) {
                  grid.rect(x, y, width, height, gp = gpar(fill = col_props(heatmap_dat[i,j]), col = NA))
                } else if (i == 7) {
                  grid.rect(x, y, width, height, gp = gpar(fill = col_age(heatmap_dat[i, j]), col = NA))
                } else {
                  grid.rect(x, y, width, height, gp = gpar(fill = col_scores(heatmap_dat[i, j]), col = NA))
                }
                grid.text(sprintf("%.1f", heatmap_dat[i, j]), x, y, gp = gpar(fontsize = 3))
              },
              row_split = c(rep("Patient", 7), 
                            "InferCNV", 
                            rep("Cytotrace", 4), 
                            rep("EMP", 3),
                            rep("HBCA", 11),
                            rep("Hallmark", 50)),
              row_gap = unit(2, "mm"),
              row_title_gp = gpar(fontsize = 5),
              row_names_gp = gpar(fontsize = 5),
              column_names_gp = gpar(fontsize = 5), 
              column_names_rot = 45,
              heatmap_legend_param = list(
                labels_gp = gpar(fontsize = 5),  # Adjust size of legend labels
                title_gp = gpar(fontsize = 5)    # Adjust size of legend title
              ),
              name = "Scores",
              border = TRUE,
              show_heatmap_legend = FALSE)
# Create custom legends
lgd_score <- Legend(col_fun = col_scores, title = "Scores")
lgd_prop <- Legend(col_fun = col_props, title = "Proportions")
lgd_age <- Legend(col_fun = col_age, title = "Age")

# Draw the heatmap with custom legends

png(file.path(PATH, "results/final_figures/epi/epi_scores_subset_clustered_heatmap_minmax_white.png"), width=1500, height=2000, res = 300)
draw(ht, heatmap_legend_list = list(lgd_score, lgd_prop, lgd_age))
dev.off()
