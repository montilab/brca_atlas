library(tidyverse)
library(randomForest)
library(ComplexHeatmap)
source(file.path(PATH, "scripts/util.R"))

PATH <- file.path(Sys.getenv("MLAB"), "projects/brcameta/brca_atlas/")
epi_data <- readRDS(file.path(PATH, "results/annotation/epi/heatmap_dat_subset.rds")) %>% t
epi_data <- epi_data %>% as.data.frame %>% rownames_to_column("Cluster")
# epi_data$k2_meta_3 <- with(epi_data, case_when(Cluster %in% c(10,4) ~ "10_4",
#                                                Cluster %in% c(1,7) ~ "1_7",
#                                                Cluster %in% c(3,8,9) ~ "3_8_9",
#                                                Cluster %in% c(2,6) ~ "2_6",
#                                                Cluster %in% c(11,12) ~ "11_12",
#                                                Cluster %in% c(0,14) ~ "0_14",
#                                                .default = Cluster))
epi_data$k2_meta_2 <- with(epi_data, case_when(Cluster %in% c(10,4,1,7) ~ "10_4_1_7",
                                               Cluster %in% c(3,8,9,2,6) ~ "3_8_9_2_6",
                                               Cluster %in% c(11,12,0,14) ~ "11_12_0_14",
                                               .default = Cluster))
# epi_data$k2_meta_1 <- with(epi_data, case_when(Cluster %in% c(10,4,1,7,5) ~ "10_4_1_7_5",
#                                                Cluster %in% c(3,8,9,2,6,13,11,12,0,14) ~ "3_8_9_2_6_13_11_12_0_14",
#                                                .default = Cluster))

# epi_data <- epi_data[1:15,]
# Assume df is your dataframe and 'group' is the column defining the groups

# Get feature importances for each group
rf_data_clusters <- list()
for(cluster in epi_data$k2_meta_2) {
  rf_epi_data <- epi_data[,-1]
  results_all <- rf_data(rf_epi_data, target_col = "k2_meta_2", target_class = cluster)
  results_imp <- results_all[["imps"]]
  results_acc <- results_all[["acc"]]
  results_imp_4 <- results_imp %>% dplyr::slice(1:4) %>% pull(feature)
  rf_data_clusters[[as.character(cluster)]][["all"]][["imp"]] <- results_imp_4
  rf_data_clusters[[as.character(cluster)]][["all"]][["direction"]] <- as.numeric(epi_data[cluster, results_imp_4]) > 0
  rf_data_clusters[[as.character(cluster)]][["all"]][["acc"]] <- results_acc
}

all_features <- colnames(epi_data)[-1]

# Create a binary matrix from the list
all_only <- lapply(rf_data_clusters, function(x) x[["all"]][["imp"]])
binary_matrix <- sapply(all_only, function(features) {
  all_features %in% features
})
rownames(binary_matrix) <- all_features
#colnames(binary_matrix) <- 0:14
heatmap_data <- binary_matrix * 1

# Barplot annotation
all_acc <- lapply(rf_data_clusters, function(x) x[["all"]][["acc"]])
all_acc <- all_acc %>% unlist %>% unname
column_ha <- HeatmapAnnotation(
  barplot = anno_barplot(all_acc),
  height = unit(6, "mm"),  
  annotation_label = "RF OOB Acc."
)

# Sorting rows
all_feats <- all_only %>% unlist %>% unname 
sorted_feats <- all_feats[!duplicated(all_feats)]
heatmap_data <- heatmap_data[sorted_feats,]
heatmap_data <- heatmap_data[rowSums(heatmap_data != 0) > 0, ]

orig_epi_data_filtered <- epi_data[,rownames(heatmap_data)] %>% t
# Create the heatmap
cell_size <- unit(3, "mm")  # Adjust this value to increase/decrease cell size

ht <- Heatmap(heatmap_data, 
        top_annotation = column_ha,
        name = "Presence", 
        border=TRUE,
        width = ncol(heatmap_data) * cell_size,
        height = nrow(heatmap_data) * cell_size,
        cell_fun = function(j, i, x, y, width, height, fill) {
          if(heatmap_data[i,j] == 0) {
            grid.rect(x = x, y = y, width = width, height = height, 
                      gp = gpar(col = "black", fill = "white"))
          } else if(orig_epi_data_filtered[i,j] < mean(orig_epi_data_filtered[i,])) {
            grid.rect(x = x, y = y, width = width, height = height, 
                      gp = gpar(col = "black", fill = "blue"))
          } else if(orig_epi_data_filtered[i,j] > mean(orig_epi_data_filtered[i,])){
            grid.rect(x = x, y = y, width = width, height = height, 
                      gp = gpar(col = "black", fill = "salmon"))
          }
        },
        row_names_gp = gpar(fontsize = 8),
        column_names_gp = gpar(fontsize = 8), 
        column_title = "Clusters", 
        row_title = "Features",
        cluster_columns = FALSE,
        cluster_rows = FALSE,
        show_row_names = TRUE, 
        show_column_names = TRUE,
        show_heatmap_legend = FALSE)
# Create the legend
lgd <- Legend(
  labels = c("High", "Low"),
  title = "Scores",
  legend_gp = gpar(fill = c("salmon", "blue")),
  grid_height = unit(6, "mm"),
  grid_width = unit(6, "mm")
)
png(file.path(PATH, "results/final_figures/epi/all_features_heatmap_k2_meta_2.png"), width=1500, height=2000, res = 300)
draw(ht, heatmap_legend_list = list(lgd))
dev.off()
