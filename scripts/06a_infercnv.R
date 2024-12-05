library(Seurat)
library(tidyverse)
library(infercnv)
library(readr)
source("util.R")

PATH <- file.path(Sys.getenv("MLAB"), "projects/brcameta/brca_atlas/")

# Loading Data
## Seurat Object
combined_seurat_rpca <- readRDS(file.path(PATH, "data/sc/combined_seurat_rpca.rds"))

batch_names <- read.delim(file.path(PATH, "scripts/06a_infercnv_data_paths"), header = FALSE)$V1
#id <- as.numeric(Sys.getenv("SGE_TASK_ID"))
id <- 7
batch_name <- batch_names[[id]]
print(batch_name)

## Genomic Position Object
gene_position <- readr::read_delim(file.path(PATH, "data/hg38_gencode_v27.txt"),
                                   delim = "\t", 
                                   escape_double = FALSE, 
                                   col_names = FALSE, 
                                   trim_ws = TRUE)
gp <- as.data.frame(gene_position)
rownames(gp) <- gp$X1
gp1 <- gp[,-1]
set.seed(123)  # For reproducibility

# Running inferCNV on subset
seurat_subset <- subset(combined_seurat_rpca, batch == batch_name)
counts_matrix <- GetAssayData(seurat_subset, slot = "counts")

rm(combined_seurat_rpca)
gc()

annotations <- data.frame(
  cell = colnames(seurat_subset),
  type = seurat_subset$cluster_broad
)
annotations$donor <- with(seurat_subset@meta.data, str_remove_all(donor, "_| "))
annotations$type_patient <- with(annotations, dplyr::case_when(type == "Epithelial" ~ paste0(type, "_", donor),
                                             .default = type))

strom_indices <- which(annotations$type == "Stromal")
imm_indices <- which(annotations$type == "Immune")
# Randomly sample 100 indices from Stromal and Immune labels
strom_to_epi_indices <- sample(strom_indices, 100)
imm_to_epi_indices <- sample(imm_indices, 100)
# Change to Non reference
strom_values <- paste0("NegControl-Strom_", annotations$donor[strom_to_epi_indices])
imm_values <- paste0("NegControl-Imm", annotations$donor[imm_to_epi_indices])
annotations$type_patient[strom_to_epi_indices] <- strom_values
annotations$type_patient[imm_to_epi_indices] <- imm_values

# Need to do this for each tumor!
# for (donor in unique(annotations$donor)) {
#   strom_indices <- which((annotations$type == "Stromal") & (annotations$donor == donor))
#   imm_indices <- which((annotations$type == "Immune") & (annotations$donor == donor))
#   # Randomly sample 100 indices from Stromal and Immune labels
#   strom_to_epi_indices <- sample(strom_indices, 100)
#   imm_to_epi_indices <- sample(imm_indices, 100)
#   # Change to Non reference
#   strom_values <- paste0("NegControl-Strom_", annotations$donor[strom_to_epi_indices])
#   imm_values <- paste0("NegControl-Imm", annotations$donor[imm_to_epi_indices])
#   annotations$type_patient[strom_to_epi_indices] <- strom_values
#   annotations$type_patient[imm_to_epi_indices] <- imm_values
# }
# print(annotations %>% 
#         dplyr::filter(str_detect(type_patient, "NegControl")) %>% 
#         dplyr::pull(type_patient) %>% 
#         table(useNA="always"))

# Save the annotations file
# annotations$cell <- NULL
annotations$type <- NULL
annotations$donor <- NULL
colnames(annotations) <- NULL
write.table(annotations, 
            file = file.path(PATH, "data/infercnv", paste0(batch_name, "_annotations.txt")), 
            sep = "\t", quote = FALSE, row.names = FALSE)

# Reference cells are pooled across tumors
# But epithelial cells are called per tumor
infercnv_obj <- infercnv::CreateInfercnvObject(raw_counts_matrix=counts_matrix,
                                    annotations_file=annotations,
                                    #delim=" ",
                                    gene_order_file=gp1,
                                    ref_group_names=c("Immune", "Stromal"))
print("Done creating infercnv obj")

out_dir <- file.path(PATH, paste0("data/infercnv/", batch_name))

if (!dir.exists(out_dir)) {
  dir.create(out_dir)
}

print("Running infercnv")
# Parameters
# 100 gene sliding window, denoised with dynamic threshold of 1.5 sd from mean
infercnv_obj <- infercnv::run(infercnv_obj,
                             cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir=out_dir,
                             cluster_by_groups=TRUE,
                             #k_obs_groups = 4,
                             denoise=TRUE,
                             num_threads = 16,
                             png_res = 100,
                             HMM=FALSE)
print(paste0("Done with infercnv for ", batch_name))
