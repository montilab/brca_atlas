library(visNetwork)
library(tidyverse)
devtools::load_all("/restricted/projectnb/brcameta/personal/andrewdr/brcasurv")
PATH <- file.path(Sys.getenv("MLAB"), "projects/brcameta/brca_atlas")
source(file.path(PATH,"../personal/andrewdr/k2_dge_script.R"))
source(file.path(PATH, "scripts/util.R"))

# Loading Data
epi_k2 <- readRDS(file.path(PATH, "results/k2/K2Taxonomer_Epi 2025_02_14_11_49_43/K2results.rds"))
strom_k2 <- readRDS(file.path(PATH, "results/k2/K2Taxonomer_Strom 2025_02_14_10_20_28/K2results.rds"))
imm_k2 <- readRDS(file.path(PATH, "results/k2/K2Taxonomer_Imm 2025_09_03_20_15_46/K2results.rds"))

# Extracting Signatures
epi_k2_dges <- getDGEresults(epi_k2)
epi_sigs <- epi_k2_dges %>%
  dplyr::mutate(list_name = paste(node, edge, direction, sep = "_")) %>%
  dplyr::ungroup() %>%
  dplyr::select(list_name, genes) %>%
  tibble::deframe()

strom_k2_dges <- getDGEresults(strom_k2)
strom_sigs <- strom_k2_dges %>%
  dplyr::mutate(list_name = paste(node, edge, direction, sep = "_")) %>%
  dplyr::ungroup() %>%
  dplyr::select(list_name, genes) %>%
  tibble::deframe()

imm_k2_dges <- getDGEresults(imm_k2)
imm_sigs <- imm_k2_dges %>%
  dplyr::mutate(list_name = paste(node, edge, direction, sep = "_")) %>%
  dplyr::ungroup() %>%
  dplyr::select(list_name, genes) %>%
  tibble::deframe()

# GSVA
# tcga_epi_gsva <- gsva_data(sigs_list = epi_sigs,
#                        adjust_prolif = TRUE,
#                        adjust_inflam = TRUE)
# metabric_epi_gsva <- gsva_data(sigs_list = epi_sigs,
#                            brca_data = "METABRIC",
#                            adjust_prolif = TRUE,
#                            adjust_inflam = TRUE)
scanb_epi_gsva <- gsva_data(sigs_list = epi_sigs,
                            brca_data = "SCANB",
                            adjust_prolif = TRUE,
                            adjust_inflam = TRUE)
# saveRDS(tcga_epi_gsva, file.path(PATH, "results/gsva/tcga_epi_gsva.rds"))
# saveRDS(metabric_epi_gsva, file.path(PATH, "results/gsva/metabric_epi_gsva.rds"))
saveRDS(scanb_epi_gsva, file.path(PATH, "results/gsva/scanb_epi_gsva.rds"))
#
# tcga_strom_gsva <- gsva_data(sigs_list = strom_sigs,
#                            adjust_prolif = TRUE,
#                            adjust_inflam = TRUE)
# metabric_strom_gsva <- gsva_data(sigs_list = strom_sigs,
#                                brca_data = "METABRIC",
#                                adjust_prolif = TRUE,
#                                adjust_inflam = TRUE)
scanb_strom_gsva <- gsva_data(sigs_list = strom_sigs,
                              brca_data = "SCANB",
                              adjust_prolif = TRUE,
                              adjust_inflam = TRUE)
# saveRDS(tcga_strom_gsva, file.path(PATH, "results/gsva/tcga_strom_gsva.rds"))
# saveRDS(metabric_strom_gsva, file.path(PATH, "results/gsva/metabric_strom_gsva.rds"))
saveRDS(scanb_strom_gsva, file.path(PATH, "results/gsva/scanb_strom_gsva.rds"))

# tcga_imm_gsva <- gsva_data(sigs_list = imm_sigs,
#                            adjust_prolif = TRUE,
#                            adjust_inflam = TRUE)
# metabric_imm_gsva <- gsva_data(sigs_list = imm_sigs,
#                                brca_data = "METABRIC",
#                                adjust_prolif = TRUE,
#                                adjust_inflam = TRUE)
scanb_imm_gsva <- gsva_data(sigs_list = imm_sigs,
                            brca_data = "SCANB",
                            adjust_prolif = TRUE,
                            adjust_inflam = TRUE)
# saveRDS(tcga_imm_gsva, file.path(PATH, "results/gsva/tcga_imm_gsva.rds"))
# saveRDS(metabric_imm_gsva, file.path(PATH, "results/gsva/metabric_imm_gsva.rds"))
saveRDS(scanb_imm_gsva, file.path(PATH, "results/gsva/scanb_imm_gsva.rds"))


