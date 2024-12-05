library(GSVA)
library(Seurat)
library(SingleCellExperiment)
library(tidyverse)
library(limma)
library("AnnotationDbi")
library("org.Hs.eg.db")
library(BiocParallel)
options(Seurat.object.assay.version = "v5")

PATH <- file.path(Sys.getenv("MLAB"), "projects/brcameta/brca_atlas/")
DATA_PATH <- file.path(Sys.getenv("CBM"), "otherStudies/scRNAseq/brca")

# Pseudobulk Data
combined_data_avg <- readRDS(file.path(PATH, "data/pseudobulk/combined_data_avg.rds"))
combined_data_avg_sce <- Seurat::as.SingleCellExperiment(combined_data_avg)

# Signature Data (IR)
# From https://pubmed.ncbi.nlm.nih.gov/20678967/
sig_data <- read.csv(file.path(PATH, "data/signatures/hardy_ir_sig.csv"), header = TRUE)
ir_up <- sig_data$UP
ir_up <- ir_up[ir_up!=""]
ir_dn <- sig_data$DOWN
ir_sigs <- list(up = ir_up, dn = ir_dn)
## Need to convert all HGNC aliases to latest first
ir_sigs <- lapply(ir_sigs, limma::alias2Symbol)

ir_gsva <- gsva(combined_data_avg_sce,
                ir_sigs,
                mx.diff=TRUE,
                verbose=FALSE)
ir_gsva$up <- assay(ir_gsva)["up",]
ir_gsva$dn <- assay(ir_gsva)["dn",]
ir_gsva$diff <- ir_gsva$up - ir_gsva$dn
saveRDS(ir_gsva, file.path(PATH, "data/pseudobulk/combined_data_hardy_gsva.rds"))

# Plotting
ir_gsva_coldata <- colData(ir_gsva) %>% as.data.frame %>% dplyr::select(c(donor,up,dn,diff))
colnames(ir_gsva_coldata) <- c("donor", "hardy_up", "hardy_dn", "hardy_diff")
fig <-  ir_gsva_coldata %>% 
          tidyr::pivot_longer(cols = !donor, names_to = "sig", values_to = "score") %>%
          ggplot() + 
          geom_violin(aes(x=sig, y=score, fill=sig)) + 
          geom_boxplot(aes(x=sig, y=score), width=0.2)
ggsave(plot=fig, filename=file.path(PATH, "results/pseudobulk/hardy_gsva_boxplot.png"), device="png")


## DE Signatures
obese_sigs <- readRDS(file.path(PATH, "brca_atlas_validation/data/sigs/obese_sigs.rds"))

# GSVA with DE signatures

obese_gsva <- gsva(combined_data_avg_sce,
                   obese_sigs,
                   mx.diff=TRUE,
                   verbose=FALSE) 
obese_gsva$all_up <- assay(obese_gsva)["all_up",]
obese_gsva$all_dn <- assay(obese_gsva)["all_dn",]
obese_gsva$all_diff <- obese_gsva$all_up - obese_gsva$all_dn
obese_gsva$epithelial_up <- assay(obese_gsva)["epi_up",]
obese_gsva$epithelial_dn <- assay(obese_gsva)["epi_dn",]
obese_gsva$epithelial_diff <- obese_gsva$epithelial_up - obese_gsva$epithelial_dn
obese_gsva$immune_up <- assay(obese_gsva)["imm_up",]
obese_gsva$immune_dn <- assay(obese_gsva)["imm_dn",]
obese_gsva$immune_diff <- obese_gsva$immune_up - obese_gsva$immune_dn
obese_gsva$endo_up <- assay(obese_gsva)["endo_up",]
obese_gsva$endo_dn <- assay(obese_gsva)["endo_dn",]
obese_gsva$endo_diff <- obese_gsva$endo_up - obese_gsva$endo_dn
obese_gsva$mesen_up <- assay(obese_gsva)["mesen_up",]
obese_gsva$mesen_dn <- assay(obese_gsva)["mesen_dn",]
obese_gsva$mesen_diff <- obese_gsva$mesen_up - obese_gsva$mesen_dn
saveRDS(obese_gsva, file.path(PATH, "data/pseudobulk/combined_data_gsva_obese.rds"))

# Storing in Metadata
ir_gsva_coldata <- colData(ir_gsva) %>% as.data.frame %>% dplyr::select(c(donor,up,dn,diff))
colnames(ir_gsva_coldata) <- c("donor", "hardy_up", "hardy_dn", "hardy_diff")
ir_gsva_coldata$donor <- str_replace_all(string=ir_gsva_coldata$donor, "-", "_")
all_metadata <- read.csv(file.path(PATH, "data/metadata/all_metadata.csv"), row.names = 1)
all_metadata <- all_metadata %>% rownames_to_column("id")
combined_metadata <- dplyr::left_join(all_metadata, ir_gsva_coldata, by = "donor")
combined_metadata <- combined_metadata %>% column_to_rownames("id")
write.csv(combined_metadata, file.path(PATH, "data/metadata/all_metadata.csv"), row.names = TRUE)

# Storing in Seurat
# combined_seurat <- readRDS(file.path(PATH, "data/sc/combined_seurat_filtered.rds"))
# ir_gsva_coldata$donor <- str_replace_all(string=ir_gsva_coldata$donor, "-", "_")
# combined_seurat@meta.data <- combined_seurat@meta.data %>% 
#   dplyr::left_join(x=., y=ir_gsva_coldata %>% as.data.frame, by="donor")
# saveRDS(combined_seurat, file.path(PATH, "data/sc/combined_seurat_filtered.rds"))
