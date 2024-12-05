library(tidyverse)
library(Seurat)
library(scCustomize)
library(parallel)

PATH <- file.path(Sys.getenv("MLAB"), "projects/brcameta/brca_atlas/")
do_save <- FALSE

# 1. Find All Markers (MAST)
imm_rpca <- readRDS(file.path(PATH, "data/sc/imm_rpca_subset.rds"))

imm_rpca$imm_broad <- with(imm_rpca@meta.data, if_else(`RNA_snn_res.0.4` %in% c("8","9","22","2",
                                                                                "3","5","11","1",
                                                                                "4","12","20","6",
                                                                                "16","19"),
                                                       "Lymphoid", 
                                                       "Myeloid"))
lymphoid_subset <- subset(imm_rpca, 
                          imm_broad == "Lymphoid")

Idents(lymphoid_subset) <- lymphoid_subset$RNA_snn_res.0.4
markers <- FindAllMarkers(lymphoid_subset,
                          assay = "RNA",
                          test.use = "MAST",
                          latent.vars = "batch",
                          only.pos = TRUE,
                          min.cells.group = 3,
                          verbose = TRUE)

# Multicore solution
# n_clust <- Idents(lymphoid_subset) %>% unique() %>% as.character() %>% as.numeric()
# mcFindMarkers <- function(i){
#   ident1 <- i
#   ident2 <- n_clust[n_clust != i]
#   table <- FindMarkers(lymphoid_subset,
#                        ident.1 = ident1, 
#                        ident.2 = ident2, 
#                        only.pos = TRUE, 
#                        test.use = "MAST", 
#                        latent.vars = "batch", 
#                        min.cells.group =3, 
#                        verbose=TRUE)
#   table$Gene.name.uniq <- rownames(table)
#   table$cluster <- rep(i, nrow(table))
#   return(table)
# }
# 
# marker_results <- list()[n_clust]
# ptm <- proc.time()
# marker_results <- parallel::mclapply(n_clust, mcFindMarkers, mc.cores = 16)
# time_diff <- proc.time() - ptm
# time_diff
# 
# # nice way to flatten list into a single DF
# dplyr::bind_rows(marker_results)
saveRDS(markers, file.path(PATH,
                           paste0("results/clustering/lymphoid_04_mast_markers.rds")))
