library(Seurat)
library(tidyverse)
PATH <- file.path(Sys.getenv("MLAB"), "projects/brcameta/brca_atlas/")
source("util.R")

wu_natgen_method <- TRUE

# Loading infercnv data
wu_natgen_2021 <- readRDS(file.path(PATH, "data/infercnv/wu_natgen_2021/run.final.infercnv_obj"))
wu_embo_2020 <- readRDS(file.path(PATH, "data/infercnv/wu_embo_2020/run.final.infercnv_obj"))
wu_genomemed_2021 <- readRDS(file.path(PATH, "data/infercnv/wu_genomemed_2021/run.final.infercnv_obj"))
pal_2021 <- readRDS(file.path(PATH, "data/infercnv/pal_2021/run.final.infercnv_obj"))
qian_2020 <- readRDS(file.path(PATH, "data/infercnv/qian_2020/run.final.infercnv_obj"))
gao_2021 <- readRDS(file.path(PATH, "data/infercnv/gao_2021/run.final.infercnv_obj"))
tietscher_2023 <- readRDS(file.path(PATH, "data/infercnv/tietscher_2023/run.final.infercnv_obj"))
bassez_2021 <- readRDS(file.path(PATH, "data/infercnv/bassez_2021/run.final.infercnv_obj"))
liu_2023 <- readRDS(file.path(PATH, "data/infercnv/liu_2023/run.final.infercnv_obj"))
wang_2024 <- readRDS(file.path(PATH, "data/infercnv/wang_2024/run.final.infercnv_obj"))

# Loading annotations 
annotation_paths <- Sys.glob(file.path(PATH, "data/infercnv/*annotations.txt"))
dataset_names <- sub(".*/(.*?)_annotations\\.txt$", 
                     "\\1", 
                     annotation_paths)
annotation_data <- lapply(annotation_paths, read.delim, header=FALSE, row.names = 1)
names(annotation_data) <- dataset_names
for(dataset in dataset_names) {
  colnames(annotation_data[[dataset]]) <- "annotation"
}

# Extracting infercnv values for observations
wu_natgen_2021mat <- wu_natgen_2021@expr.data[, unlist(wu_natgen_2021@observation_grouped_cell_indices)]
wu_embo_2020mat <- wu_embo_2020@expr.data[, unlist(wu_embo_2020@observation_grouped_cell_indices)]
wu_genomemed_2021mat <- wu_genomemed_2021@expr.data[, unlist(wu_genomemed_2021@observation_grouped_cell_indices)]
pal_2021mat <- pal_2021@expr.data[, unlist(pal_2021@observation_grouped_cell_indices)]
qian_2020mat <- qian_2020@expr.data[, unlist(qian_2020@observation_grouped_cell_indices)]
gao_2021mat <- gao_2021@expr.data[, unlist(gao_2021@observation_grouped_cell_indices)]
tietscher_2023mat <- tietscher_2023@expr.data[, unlist(tietscher_2023@observation_grouped_cell_indices)]
bassez_2021mat <- bassez_2021@expr.data[, unlist(bassez_2021@observation_grouped_cell_indices)]
liu_2023mat <- liu_2023@expr.data[, unlist(liu_2023@observation_grouped_cell_indices)]
wang_2024mat <- wang_2024@expr.data[, unlist(wang_2024@observation_grouped_cell_indices)]

if(wu_natgen_method) {
  # Mean of squared, scaled cnv scores
  
  # Scaling each gene to -1,1
  wu_natgen_2021mat_scaled <- t(scale_to_range(t(wu_natgen_2021mat)))
  wu_embo_2020mat_scaled <- t(scale_to_range(t(wu_embo_2020mat)))
  wu_genomemed_2021mat_scaled <- t(scale_to_range(t(wu_genomemed_2021mat)))
  pal_2021mat_scaled <- t(scale_to_range(t(pal_2021mat)))
  qian_2020mat_scaled <- t(scale_to_range(t(qian_2020mat)))
  gao_2021mat_scaled <- t(scale_to_range(t(gao_2021mat)))
  tietscher_2023mat_scaled <- t(scale_to_range(t(tietscher_2023mat)))
  bassez_2021mat_scaled <- t(scale_to_range(t(bassez_2021mat)))
  liu_2023mat_scaled <- t(scale_to_range(t(liu_2023mat)))
  wang_2024mat_scaled <- t(scale_to_range(t(wang_2024mat)))
  
  # Squaring and taking mean
  wu_natgen_2021sum <- colSums(wu_natgen_2021mat_scaled^2) / nrow(wu_natgen_2021mat_scaled)
  wu_embo_2020sum <- colSums(wu_embo_2020mat_scaled^2) / nrow(wu_embo_2020mat_scaled)
  wu_genomemed_2021sum <- colSums(wu_genomemed_2021mat_scaled^2) / nrow(wu_genomemed_2021mat_scaled)
  pal_2021sum <- colSums(pal_2021mat_scaled^2) / nrow(pal_2021mat_scaled)
  qian_2020sum <- colSums(qian_2020mat_scaled^2) / nrow(qian_2020mat_scaled)
  gao_2021sum <- colSums(gao_2021mat_scaled^2) / nrow(gao_2021mat_scaled)
  tietscher_2023sum <- colSums(abs(tietscher_2023mat_scaled)) / nrow(tietscher_2023mat_scaled)
  bassez_2021sum <- colSums(bassez_2021mat_scaled^2) / nrow(bassez_2021mat_scaled)
  liu_2023sum <- colSums(liu_2023mat_scaled^2) / nrow(liu_2023mat_scaled)
  wang_2024sum <- colSums(wang_2024mat_scaled^2) / nrow(wang_2024mat_scaled)  
} else {
  
  # Just taking mean of absolute value of scores
  wu_natgen_2021sum <- colSums(abs(wu_natgen_2021mat)) / nrow(wu_natgen_2021mat)
  wu_embo_2020sum <- colSums(abs(wu_embo_2020mat)) / nrow(wu_embo_2020mat)
  wu_genomemed_2021sum <- colSums(abs(wu_genomemed_2021mat)) / nrow(wu_genomemed_2021mat)
  pal_2021sum <- colSums(abs(pal_2021mat)) / nrow(pal_2021mat)
  qian_2020sum <- colSums(abs(qian_2020mat)) / nrow(qian_2020mat)
  gao_2021sum <- colSums(abs(gao_2021mat)) / nrow(gao_2021mat)
  tietscher_2023sum <- colSums(abs(tietscher_2023mat)) / nrow(tietscher_2023mat)
  bassez_2021sum <- colSums(abs(bassez_2021mat)) / nrow(bassez_2021mat)
  liu_2023sum <- colSums(abs(liu_2023mat)) / nrow(liu_2023mat)
  wang_2024sum <- colSums(abs(wang_2024mat)) / nrow(wang_2024mat)
}

# par(mfrow=c(5,2))
# hist(wu_natgen_2021sum, breaks=50, freq=FALSE, main = "wu_natgen_2021")
# hist(wu_embo_2020sum, breaks=50, freq=FALSE, main = "wu_embo_2020")
# hist(wu_genomemed_2021sum, breaks=50, freq=FALSE, main = "wu_genomemed_2021")
# hist(pal_2021sum, breaks=50, freq=FALSE, main = "pal_2021")
# hist(qian_2020sum, breaks=50, freq=FALSE, main = "qian_2020")
# hist(gao_2021sum, breaks=50, freq=FALSE, main = "gao_2021")
# hist(tietscher_2023sum, breaks=50, freq=FALSE, main = "tietscher_2023")
# hist(bassez_2021sum, breaks=50, freq=FALSE, main = "bassez_2021")
# hist(liu_2023sum, breaks=50, freq=FALSE, main = "liu_2023")
# hist(wang_2024sum, breaks=50, freq=FALSE, main = "wang_2024")

# Plot malignancy scores for each dataset along with controls
# Epithelial-patient, neg control- patient
wu_natgen_2021_malig_data <- wu_natgen_2021sum %>% as.data.frame()
colnames(wu_natgen_2021_malig_data) <- "score"
wu_natgen_2021_malig_data <- merge(wu_natgen_2021_malig_data, annotation_data[["wu_natgen_2021"]], by=0)
colnames(wu_natgen_2021_malig_data) <- c("cell_id", "score", "annotation")
wu_natgen_2021_malig_data$batch <- "wu_natgen_2021"

wu_embo_2020_malig_data <- wu_embo_2020sum %>% as.data.frame()
colnames(wu_embo_2020_malig_data) <- "score"
wu_embo_2020_malig_data <- merge(wu_embo_2020_malig_data, annotation_data[["wu_embo_2020"]], by=0)
colnames(wu_embo_2020_malig_data) <- c("cell_id", "score", "annotation")
wu_embo_2020_malig_data$batch <- "wu_embo_2020"

wu_genomemed_2021_malig_data <- wu_genomemed_2021sum %>% as.data.frame()
colnames(wu_genomemed_2021_malig_data) <- "score"
wu_genomemed_2021_malig_data <- merge(wu_genomemed_2021_malig_data, annotation_data[["wu_genomemed_2021"]], by=0)
colnames(wu_genomemed_2021_malig_data) <- c("cell_id", "score", "annotation")
wu_genomemed_2021_malig_data$batch <- "wu_genomemed_2021"

pal_2021_malig_data <- pal_2021sum %>% as.data.frame()
colnames(pal_2021_malig_data) <- "score"
pal_2021_malig_data <- merge(pal_2021_malig_data, annotation_data[["pal_2021"]], by=0)
colnames(pal_2021_malig_data) <- c("cell_id", "score", "annotation")
pal_2021_malig_data$batch <- "pal_2021"

qian_2020_malig_data <- qian_2020sum %>% as.data.frame()
colnames(qian_2020_malig_data) <- "score"
qian_2020_malig_data <- merge(qian_2020_malig_data, annotation_data[["qian_2020"]], by=0)
colnames(qian_2020_malig_data) <- c("cell_id", "score", "annotation")
qian_2020_malig_data$batch <- "qian_2020"

gao_2021_malig_data <- gao_2021sum %>% as.data.frame()
colnames(gao_2021_malig_data) <- "score"
gao_2021_malig_data <- merge(gao_2021_malig_data, annotation_data[["gao_2021"]], by=0)
colnames(gao_2021_malig_data) <- c("cell_id", "score", "annotation")
gao_2021_malig_data$batch <- "gao_2021"

tietscher_2023_malig_data <- tietscher_2023sum %>% as.data.frame()
colnames(tietscher_2023_malig_data) <- "score"
tietscher_2023_malig_data <- merge(tietscher_2023_malig_data, annotation_data[["tietscher_2023"]], by=0)
colnames(tietscher_2023_malig_data) <- c("cell_id", "score", "annotation")
tietscher_2023_malig_data$batch <- "tietscher_2023"

bassez_2021_malig_data <- bassez_2021sum %>% as.data.frame()
colnames(bassez_2021_malig_data) <- "score"
bassez_2021_malig_data <- merge(bassez_2021_malig_data, annotation_data[["bassez_2021"]], by=0)
colnames(bassez_2021_malig_data) <- c("cell_id", "score", "annotation")
bassez_2021_malig_data$batch <- "bassez_2021"

liu_2023_malig_data <- liu_2023sum %>% as.data.frame()
colnames(liu_2023_malig_data) <- "score"
liu_2023_malig_data <- merge(liu_2023_malig_data, annotation_data[["liu_2023"]], by=0)
colnames(liu_2023_malig_data) <- c("cell_id", "score", "annotation")
liu_2023_malig_data$batch <- "liu_2023"

wang_2024_malig_data <- wang_2024sum %>% as.data.frame()
colnames(wang_2024_malig_data) <- "score"
wang_2024_malig_data <- merge(wang_2024_malig_data, annotation_data[["wang_2024"]], by=0)
colnames(wang_2024_malig_data) <- c("cell_id", "score", "annotation")
wang_2024_malig_data$batch <- "wang_2024"
combined_malig_data <- rbind(
  wu_natgen_2021_malig_data,
  wu_embo_2020_malig_data,
  wu_genomemed_2021_malig_data,
  pal_2021_malig_data,
  qian_2020_malig_data,
  gao_2021_malig_data,
  tietscher_2023_malig_data,
  bassez_2021_malig_data,
  liu_2023_malig_data,
  wang_2024_malig_data
)

infercnv_plots <- list()
for(i in seq_along(dataset_names)) {
  dataset_name <- dataset_names[[i]]
  batch_data <- combined_malig_data %>% dplyr::filter(batch == dataset_name)
  y_95 <- unique(batch_data$norm_95)
  y_5 <- unique(batch_data$norm_5)
  p <- batch_data %>%
    dplyr::mutate(cancer = str_detect(annotation, "Epithelial")) %>%
    ggplot() + 
    geom_boxplot(aes(x = annotation, y = score, fill = cancer)) + 
    geom_hline(color = "maroon", linetype = "dashed", yintercept = y_95) +
    geom_hline(color = "maroon", linetype = "dashed", yintercept = y_5) +
    scale_fill_manual(values = c("TRUE"="lightsalmon2", "FALSE"="lightgoldenrod")) +
    labs(title = dataset_name,
         y = "Malignancy Score",
         x = "Sample") +
    theme(legend.position = "right",
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16))
  infercnv_plots[[i]] <- p
  ggsave(filename = file.path(PATH, paste0("results/infercnv/", dataset_name, "_scores.png")), width = 7, height = 5)  
}

library(cowplot)
plot_grid(plotlist = infercnv_plots, ncol = 2)
ggsave(filename = file.path(PATH, "results/infercnv/all_data_scores.png"), width = 11, height = 14)  
# Add data to seurat object
# Plot malignancy scores for each epithelial cluster
#epi_rpca <- readRDS(file.path(PATH, "data/sc/epi_rpca_subset.rds"))
rownames(combined_malig_data) <- combined_malig_data$cell_id
saveRDS(combined_malig_data, file.path(PATH, "data/infercnv/combined_infercnv_data.rds"))

combined_malig_data <- merge(combined_malig_data, epi_rpca@meta.data, all.x = TRUE, by = 0)
rownames(combined_malig_data) <- combined_malig_data$Row.names
combined_malig_data$Row.names <- NULL
combined_malig_data$score %>% fivenum(na.rm=TRUE)
combined_malig_data$group <- with(combined_malig_data, case_when(str_detect(annotation, "Epithelial") ~ `RNA_snn_res.0.2`,
                                                            str_detect(annotation, "Strom") ~ "Stromal",
                                                            str_detect(annotation, "Imm") ~ "Immune"))
combined_malig_data %>% 
  dplyr::mutate(cancer = str_detect(annotation, "Epithelial")) %>%
  ggplot() + 
  geom_boxplot(aes(x = reorder(group, -score, mean), y = score, fill = cancer)) + 
  scale_fill_manual(values = c("TRUE"="lightsalmon2", "FALSE"="lightgoldenrod")) +
  labs(title = "Malignancy Scores by Cluster",
       y = "Malignancy Score",
       x = "Sample") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        axis.ticks.x = element_blank(),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16))
ggsave(filename = file.path(PATH, paste0("results/infercnv/all_clust_squared_scaled_scores.png")), width = 7, height = 8)  

# Plot malignancy scores for each subtype
combined_malig_data$group <- with(combined_malig_data, case_when(str_detect(annotation, "Epithelial") ~ subtype_new,
                                                                 str_detect(annotation, "Strom") ~ "Stromal",
                                                                 str_detect(annotation, "Imm") ~ "Immune"))
combined_malig_data %>% 
  dplyr::mutate(cancer = str_detect(annotation, "Epithelial")) %>%
  ggplot() + 
  geom_boxplot(aes(x = reorder(group, -score, mean), y = score, fill = cancer)) + 
  scale_fill_manual(values = c("TRUE"="lightsalmon2", "FALSE"="lightgoldenrod")) +
  labs(title = "Malignancy Scores by Subtype",
       y = "Malignancy Score",
       x = "Sample") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        axis.ticks.x = element_blank(),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16))
ggsave(filename = file.path(PATH, paste0("results/infercnv/subtype_squared_scaled_scores.png")), width = 7, height = 8)  

# Saving to epi_rpca

saveRDS(combined_malig_data, file.path(PATH, "data/infercnv/combined_infercnv_data.rds"))

# Pal 2021 experiment (sample random 100 epi cells)
dataset_name == "pal_2021"
sampled_data <- combined_malig_data %>%
    dplyr::filter((batch.x == dataset_name) &
                     str_detect(annotation, "Epithelial")) %>%
    dplyr::group_by(annotation) %>% 
    dplyr::slice_sample(n = 100)
sampled_data$annotation <- str_replace(sampled_data$annotation, "Epithelial", "Sampled")
combined_malig_data <- rbind(combined_malig_data, sampled_data)

combined_malig_data %>% 
  dplyr::filter(batch.x == dataset_name) %>%
  dplyr::filter(str_detect(annotation, "Epithelial|Sampled")) %>%
  dplyr::mutate(cancer = str_detect(annotation, "Epithelial")) %>%
  ggplot() + 
  geom_boxplot(aes(x = annotation, y = score, fill = cancer)) + 
  scale_fill_manual(values = c("TRUE"="lightsalmon2", "FALSE"="lightgoldenrod")) +
  labs(title = dataset_name,
       y = "Malignancy Score",
       x = "Sample") +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16))
ggsave(filename = file.path(PATH, paste0("results/infercnv/", dataset_name, "_sampled_squared_scores.png")), width = 7, height = 5) 


# Finding the amount of cells per tumor outside the normal envelope

combined_malig_data <- readRDS(file.path(PATH, "data/infercnv/combined_infercnv_data.rds"))
normal_malig_ranges <- combined_malig_data %>%
  dplyr::group_by(batch) %>% 
  dplyr::filter(str_detect(annotation, "NegControl")) %>%
  dplyr::summarise(
    norm_25 = quantile(score, 0.025),
    norm_5 = quantile(score, 0.05),
    norm_95 = quantile(score, 0.95),
    norm_975 = quantile(score, 0.975)
  )
combined_malig_data <- combined_malig_data %>% 
  dplyr::left_join(x=., y=normal_malig_ranges, by="batch") %>%
  dplyr::mutate(outside_norm_95 = case_when(score > norm_975 ~ TRUE,
                                         score < norm_25 ~ TRUE,
                                         .default = FALSE),
                outside_norm_90 = case_when(score > norm_95 ~ TRUE,
                                            score < norm_5 ~ TRUE,
                                            .default = FALSE))
saveRDS(combined_malig_data, file.path(PATH, "data/infercnv/combined_infercnv_data.rds"))

