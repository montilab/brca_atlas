library(Seurat)
PATH <- file.path(Sys.getenv("MLAB"), "projects/brcameta/brca_atlas")

# Loading Data
combined_seurat <- readRDS(file=file.path(PATH, "data/sc/combined_seurat_filtered.rds"))

# Averaging by donor (sample)
# combined_data_avg <- Seurat::AggregateExpression(combined_seurat, group.by="donor", return.seurat=TRUE)
combined_data_avg <- Seurat::AverageExpression(combined_seurat, group.by="donor", return.seurat=TRUE)
# Finding differing amounts of variable features for model training later
n_features <- c(2000,5000,10000)
for(n_feature in n_features) {
  combined_seurat_avg <- Seurat::FindVariableFeatures(combined_data_avg, nfeatures=n_feature)
  var_genes <- Seurat::VariableFeatures(combined_seurat_avg)
  saveRDS(var_genes, file.path(PATH, paste0("/data/var_genes/combined_var_avg_", n_feature, ".rds")))
}

saveRDS(combined_data_avg, file.path(PATH, "data/pseudobulk/combined_data_avg.rds"))