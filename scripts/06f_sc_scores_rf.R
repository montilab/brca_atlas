library(tidyverse)
library(randomForest)
PATH <- file.path(Sys.getenv("MLAB"), "projects/brcameta/brca_atlas/")
source(file.path(PATH, "scripts/util.R"))
library(doParallel)
library(foreach)
library(fastDummies)
# Set up parallel backend
# num_cores <- detectCores() - 1  # Use all cores except one
registerDoParallel(cores = 8)

epi_data <- readRDS(file.path(PATH, "results/annotation/epi/scores_cell.rds"))
epi_data <- epi_data %>% dplyr::filter(Cluster %in% 0:14)
epi_data$Cluster <- droplevels(epi_data$Cluster)
epi_data$Cytotrace_Class <- droplevels(epi_data$Cytotrace_Class)
epi_data$Subtype <- NULL
epi_data$`Cell id` <- NULL
epi_data$InferCNV_Malig <- epi_data$InferCNV_Malig*1
epi_data <- dummy_cols(epi_data, select_columns = c("Pam50", "Cytotrace_Class"), 
                       remove_selected_columns = TRUE)
colnames(epi_data) <- str_replace_all(colnames(epi_data), "Cytotrace_Class_", replacement = "Cyto. ")
colnames(epi_data) <- str_replace_all(colnames(epi_data), "Cytotrace_Score", replacement = "Cyto. Score")
colnames(epi_data) <- str_replace_all(colnames(epi_data), "Pam50_", replacement = "Pam50 ")


# Random Forest needs complete data
epi_data <- epi_data[complete.cases(epi_data),]
table(epi_data$Cluster)

# 1. Random Forest for all clusters
rf_data_clusters <- foreach(cluster = unique(epi_data$Cluster), .combine = c) %dopar% {
  results_all <- rf_data(epi_data, target_col = "Cluster", target_class = cluster)
  results_imp <- results_all[["imps"]]
  results_acc <- results_all[["acc"]]
  results_imp_feats <- results_imp %>% dplyr::slice(1:5) %>% pull(feature)
  
  cluster_result <- list()
  cluster_result[[as.character(cluster)]] <- list(
    all = list(
      imp = results_imp_feats,
      acc = results_acc
    )
  )
  
  print(paste0("Done with ", cluster))
  cluster_result
}

saveRDS(rf_data_clusters, file.path(PATH, "results/annotation/epi/rf_all_clusters.rds"))

# 1. Random Forest for level 3 clusters
epi_data$k2_meta_3 <- with(epi_data, case_when(Cluster %in% c(10,4) ~ "10_4",
                                               Cluster %in% c(1,7) ~ "1_7",
                                               Cluster %in% c(3,8,9) ~ "3_8_9",
                                               Cluster %in% c(2,6) ~ "2_6",
                                               Cluster %in% c(11,12) ~ "11_12",
                                               Cluster %in% c(0,14) ~ "0_14",
                                               .default = Cluster))
table(epi_data$k2_meta_3)

# Get feature importances for each group
rf_data_clusters <- foreach(cluster = unique(epi_data$k2_meta_3), .combine = c) %dopar% {
  results_all <- rf_data(epi_data %>% dplyr::select(-Cluster),
                         target_col = "k2_meta_3",
                         target_class = cluster)
  results_imp <- results_all[["imps"]]
  results_acc <- results_all[["acc"]]
  results_imp_feats <- results_imp %>% dplyr::slice(1:5) %>% pull(feature)

  cluster_result <- list()
  cluster_result[[as.character(cluster)]] <- list(
    all = list(
      imp = results_imp_feats,
      acc = results_acc
    )
  )

  print(paste0("Done with ", cluster))
  cluster_result
}

saveRDS(rf_data_clusters, file.path(PATH, "results/annotation/epi/rf_k2_3_clusters.rds"))

# 1. Random Forest for level 2 clusters
epi_data <- epi_data %>% dplyr::select(-k2_meta_3)
epi_data$k2_meta_2 <- with(epi_data, case_when(Cluster %in% c(10,4,1,7) ~ "10_4_1_7",
                                               Cluster %in% c(3,8,9,2,6) ~ "3_8_9_2_6",
                                               Cluster %in% c(11,12,0,14) ~ "11_12_0_14",
                                               .default = Cluster))

table(epi_data$k2_meta_2)

# Get feature importances for each group
rf_data_clusters <- foreach(cluster = unique(epi_data$k2_meta_2), .combine = c) %dopar% {
  results_all <- rf_data(epi_data %>% dplyr::select(-Cluster),
                         target_col = "k2_meta_2",
                         target_class = cluster)
  results_imp <- results_all[["imps"]]
  results_acc <- results_all[["acc"]]
  results_imp_feats <- results_imp %>% dplyr::slice(1:5) %>% pull(feature)

  cluster_result <- list()
  cluster_result[[as.character(cluster)]] <- list(
    all = list(
      imp = results_imp_feats,
      acc = results_acc
    )
  )

  print(paste0("Done with ", cluster))
  cluster_result
}

saveRDS(rf_data_clusters, file.path(PATH, "results/annotation/epi/rf_k2_2_clusters.rds"))

# 4. Random Forest Level 1 clusters
epi_data <- epi_data %>% dplyr::select(-k2_meta_2)
epi_data$k2_meta_1 <- with(epi_data, case_when(Cluster %in% c(10,4,1,7,5) ~ "10_4_1_7_5",
                                               Cluster %in% c(3,8,9,2,6,13,11,12,0,14) ~ "3_8_9_2_6_13_11_12_0_14",
                                               .default = Cluster))
table(epi_data$k2_meta_1)

# Get feature importances for each group
rf_data_clusters <- foreach(cluster = unique(epi_data$k2_meta_1), .combine = c) %dopar% {
  results_all <- rf_data(epi_data %>% dplyr::select(-Cluster),
                         target_col = "k2_meta_1",
                         target_class = cluster)
  results_imp <- results_all[["imps"]]
  results_acc <- results_all[["acc"]]
  results_imp_feats <- results_imp %>% dplyr::slice(1:5) %>% pull(feature)

  cluster_result <- list()
  cluster_result[[as.character(cluster)]] <- list(
    all = list(
      imp = results_imp_feats,
      acc = results_acc
    )
  )

  print(paste0("Done with ", cluster))
  cluster_result
}

saveRDS(rf_data_clusters, file.path(PATH, "results/annotation/epi/rf_k2_1_clusters.rds"))
