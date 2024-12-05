library(stringr)

singler_short <- function(singler_annot) {
  # Epithelial
  singler_annot[str_detect(singler_annot, "luminal epithelial")] <- "Epithelial"
  singler_annot[str_detect(singler_annot, "mammary gland epithelial cell")] <- "Epithelial"
  singler_annot[str_detect(singler_annot, "basal cell")] <- "Epithelial"
  # Immune
  singler_annot[str_detect(singler_annot, "CD4")] <- "CD4 T-Cells"
  singler_annot[str_detect(singler_annot, "CD8")] <- "CD8 T-Cells"
  singler_annot[str_detect(singler_annot, "gamma-delta")] <- "Gamma-Delta"
  singler_annot[str_detect(singler_annot, "T cell")] <- "T-Cells"
  singler_annot[str_detect(singler_annot, "regulatory")] <- "T-Regs"
  singler_annot[str_detect(singler_annot, "natural killer") | str_detect(singler_annot, "NK")] <- "NK Cells "
  singler_annot[str_detect(singler_annot, "memory B cell")] <- "B Cells"
  singler_annot[str_detect(singler_annot, "naive B cell")] <- "B Cells"
  singler_annot[str_detect(singler_annot, "plasma")] <- "Plasma"
  singler_annot[str_detect(singler_annot, "macrophage")] <- "Macrophages"
  singler_annot[str_detect(singler_annot, "monocyte")] <- "Monocytes"
  singler_annot[str_detect(singler_annot, "dendritic")] <- "Dendritic"
  # Stromal
  singler_annot[str_detect(singler_annot, "endothelial") | str_detect(singler_annot, "Endothelial")] <- "Endothelial"
  singler_annot[str_detect(singler_annot, "capillary")] <- "Endothelial"
  singler_annot[str_detect(singler_annot, "muscle")] <- "Smooth Muscle"
  
  singler_annot <- str_to_title(singler_annot)
  return(singler_annot)
}

singler_broad <- function(singler_annot) {
  # Epi
  singler_annot[str_detect(singler_annot, "epithelial|basal")] <- "Epithelial" 
  # Imm
  singler_annot[str_detect(singler_annot, 
                           "T cell|B cell|myeloid|macrophage|neutrophil|killer|lymphocyte|monocyte|plasma|dendritic|mast")] <- "Immune"
  # Strom
  singler_annot[str_detect(singler_annot, "endo|capillary|muscle|fibro|pericyte")] <- "Stromal"
  return(singler_annot)
}

celltypist_short <- function(celltypist_annot) {
  # Epithelial
  celltypist_annot[str_detect(celltypist_annot, "basal")] <- "Epithelial"
  celltypist_annot[str_detect(celltypist_annot, "LummHR")] <- "Epithelial"
  celltypist_annot[str_detect(celltypist_annot, "Lumsec")] <- "Epithelial"
  # Immune
  celltypist_annot[str_detect(celltypist_annot, "CD4")] <- "CD4 T-Cells"
  celltypist_annot[str_detect(celltypist_annot, "CD8")] <- "CD8 T-Cells"
  celltypist_annot[str_detect(celltypist_annot, "DC")] <- "Dendritic"
  celltypist_annot[str_detect(celltypist_annot, "GD")] <- "Gamma-Delta"
  celltypist_annot[str_detect(celltypist_annot, "NK")] <- "NK Cells "
  celltypist_annot[str_detect(celltypist_annot, "bmem")] <- "B Cells"
  celltypist_annot[str_detect(celltypist_annot, "b_naive")] <- "B Cells"
  celltypist_annot[str_detect(celltypist_annot, "plasma")] <- "Plasma"
  celltypist_annot[str_detect(celltypist_annot, "Macro")] <- "Macrophages"
  celltypist_annot[str_detect(celltypist_annot, "Mono")] <- "Monocytes"
  celltypist_annot[str_detect(celltypist_annot, "mye-prol")] <- "Myeloid Prolif."
  celltypist_annot[str_detect(celltypist_annot, "T_prol")] <- "Lymphoid Prolif."
  # Stromal
  celltypist_annot[str_detect(celltypist_annot, "Fibro")] <- "Fibroblasts"
  celltypist_annot[str_detect(celltypist_annot, "pericytes")] <- "Pericytes"
  celltypist_annot[str_detect(celltypist_annot, "Vas")] <- "Endothelial"
  celltypist_annot[str_detect(celltypist_annot, "vsmc")] <- "Smooth Muscle"
  
  
  celltypist_annot <- str_to_title(celltypist_annot)
  return(celltypist_annot)
}

celltypist_broad <- function(celltypist_annot) {
  # Epi
  celltypist_annot[str_detect(celltypist_annot, "LummHR|basal|Lumsec")] <- "Epithelial" 
  # Imm
  celltypist_annot[str_detect(celltypist_annot, "CD4|CD8|DC|GD|NK|bmem|b_naive|plasma|Macro|Mono|mye-prol|T_prol|Mast|Neutrophil")] <- "Immune"
  # Strom
  celltypist_annot[str_detect(celltypist_annot, "Fibro|pericytes|Vas|vsmc|Lymph")] <- "Stromal"
  return(celltypist_annot)
}

author_short <- function(author_annot) {
  # Epithelial
  author_annot[str_detect(author_annot, "Cancer")] <- "Epithelial"
  author_annot[str_detect(author_annot, "epithelial")] <- "Epithelial"
  # Immune
  author_annot[str_detect(author_annot, "CD4")] <- "CD4 T-Cells"
  author_annot[str_detect(author_annot, "CD8")] <- "CD8 T-Cells"
  author_annot[str_detect(author_annot, "T-cells")] <- "T-Cells"
  author_annot[str_detect(author_annot, "Tfh")] <- "T-Cells"
  author_annot[str_detect(author_annot, "DC")] <- "Dendritic"
  author_annot[str_detect(author_annot, "NK")] <- "Nk Cells"
  author_annot[str_detect(author_annot, "blasts")] <- "Plasma"
  author_annot[str_detect(author_annot, "granu")] <- "Neutrophil"
  author_annot[str_detect(author_annot, "Monocytes")] <- "Monocytes"
  author_annot[str_detect(author_annot, "Macrophages Cycling")] <- "Myeloid prol."
  author_annot[str_detect(author_annot, "T-Cells Cycling")] <- "Lymphoid prol."
  # Stromal
  author_annot[str_detect(author_annot, "PVL")] <- "Pericyte"
  author_annot[str_detect(author_annot, "CAFs")] <- "CAFs"
  author_annot[str_detect(author_annot, "EC")] <- "Endothelial"
  
  author_annot <- str_to_title(author_annot)
  return(author_annot)
}

author_new <- function(author_annot) {
  author_annot[is.na(author_annot)] <- "Unassigned"
  #Stromal
  author_annot[author_annot == "endothelial"] <- "Endothelial"
  author_annot[author_annot == "Endothelial_cell"] <- "Endothelial"
  author_annot[author_annot == "PVL cells"] <- "PVL"
  author_annot[author_annot == "fibroblast"] <- "Fibroblast"
  #Immune
  author_annot[author_annot == "B cells"] <- "B-cells"
  author_annot[author_annot == "Plasma_Cells"] <- "Plasma"
  author_annot[author_annot == "B_Cells"] <- "B-cells"
  author_annot[author_annot == "B_cell"] <- "B-cells"
  author_annot[author_annot == "B cell"] <- "B-cells"
  author_annot[author_annot == "myeloid"] <- "Myeloid"
  author_annot[author_annot == "Mast_cell"] <- "Mast"
  author_annot[author_annot == "Myeloid_cell"] <- "Myeloid"
  author_annot[author_annot == "pDCs"] <- "pDC"
  author_annot[author_annot == "plasma cell"] <- "Plasma"
  author_annot[author_annot == "Plasma_Cells"] <- "Plasma"
  author_annot[author_annot == "T_cell" | author_annot == "T_cells_unassigned"] <- "T-cells"
  #Epi
  author_annot[author_annot == "Normal Epithelial"] <- "Epithelial"
  author_annot[author_annot == "epithelial"] <- "Epithelial"
  author_annot[author_annot == "Malignant"] <- "Epithelial"
  author_annot[author_annot == "Cancer_cell"] <- "Epithelial"
  author_annot[author_annot == "Cancer/Epithelial"] <- "Epithelial"
  author_annot[author_annot == "Cancer/Epithelial Cycling"] <- "Epithelial"
  author_annot[author_annot == "Cancer Epithelial"] <- "Epithelial"
  author_annot[author_annot == "Cancer Epithelial Cycling"] <- "Epithelial"
  author_annot[author_annot == "Epithelial_Basal"] <- "Epithelial"
  author_annot[author_annot == "Epithelial_Basal_Cycling"] <- "Epithelial"
  author_annot[author_annot == "Epithelial_Luminal_Mature"] <- "Epithelial"
  
  return(author_annot)
}

author_broad <- function(author_annot) {
  # Epi
  author_annot[str_detect(author_annot, "Cancer|Epithelial|epithelial|Malignant")] <- "Epithelial" 
  # Imm
  author_annot[str_detect(author_annot, 
                          "DC|T-cell|granu|Lymph|B-cell|NK|Mono|Macro|Tfh|Plasma|Mast|Dendritic|T-Regs|Myeloid")] <- "Immune"
  # Strom
  author_annot[str_detect(author_annot, "Pericyte|EC|Endo|Fibro|PVL|CAF")] <- "Stromal"
  return(author_annot)
}

dotplot_ident <- function(annot) {
  # Epithelial
  # Immune
  annot[str_detect(annot, "Gamma-Delta")] <- "Lymphoid Prol."
  # Stromal
  annot[str_detect(annot, "EC|Lymph|Pericytes")] <- "Endothelial"
  
  annot <- str_to_title(annot)
  return(annot)
}

scale_to_range <- function(x, new_min = -1, new_max = 1) {
  x_min <- apply(x, 2, min)
  x_max <- apply(x, 2, max)
  t((t(x) - x_min) / (x_max - x_min) * (new_max - new_min) + new_min)
}


sccomp_prep_data <- function(sccomp_df, param, v_c = "c", fdr = 0.05, order = NULL) {
  
  stopifnot(all(c("parameter", "c_FDR", "c_effect", "cluster_annot") %in% colnames(sccomp_df)))
  
  if (v_c == "c") {
    sccomp_df <- sccomp_df %>% 
      dplyr::filter(parameter == param) %>%
      dplyr::mutate(diffexpressed = case_when((c_FDR < fdr) & (c_effect < 0) ~ "dn",
                                              (c_FDR < fdr) & (c_effect > 0) ~ "up",
                                              .default = "ns"),
                    y = rank(c_effect))
  } else if (v_c == "v") {
    sccomp_df <- sccomp_df %>% 
      dplyr::filter(parameter == param) %>%
      dplyr::mutate(diffexpressed = case_when((v_FDR < fdr) & (v_effect < 0) ~ "dn",
                                              (v_FDR < fdr) & (v_effect > 0) ~ "up",
                                              .default = "ns"),
                    y = rank(v_effect))
  }

  if (!is.null(order)) {
    sccomp_df$cluster_annot <- factor(sccomp_df$cluster_annot, levels = order)
    sccomp_df <- sccomp_df[order(sccomp_df$cluster_annot),]
  }
  return(sccomp_df)
}

sccomp_plot_data <- function(sccomp_df,
                             v_c = "c",
                             up = "#FF9999",
                             dn = "#56B4E9",
                             sort_y = TRUE) {
  
  stopifnot(all(c("c_FDR", "c_effect", "diffexpressed", "c_upper", "c_lower", "cluster_annot") %in% colnames(sccomp_df)))
  
  if(!sort_y) {
    sccomp_df$y <- nrow(sccomp_df):1
  }
  pd <- ggplot2::position_dodge(0.1)
  
  if (v_c == "c") {
    p <- sccomp_df %>%
      ggplot2::ggplot(aes(x = c_effect, y = y, col = diffexpressed, size=c_FDR)) +
      ggplot2::geom_vline(xintercept = 0, col = "gray") +
      ggplot2::geom_point() +
      ggplot2::geom_errorbar(
        aes(xmin = c_lower, xmax = c_upper),
        width = 0.5, position = pd, orientation = "y", size=0.5
      ) +
      ggplot2::scale_y_continuous(
        breaks = sccomp_df |> dplyr::pull(y),
        labels = sccomp_df |> dplyr::pull("cluster_annot")) +
      ggplot2::scale_size_continuous(trans = ggforce::trans_reverser("log10")) +
      ggplot2::xlab("Credible Interval") +
      ggplot2::theme_classic()  +
      ggplot2::scale_color_manual(values = c("dn" = dn, "ns" = "gray", "up" = up)) +
      ggplot2::theme(
        axis.title.x = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        axis.title.y = element_blank() ,
        axis.text.y = element_text(size = 14),
      ) +
      ggplot2::labs(color = "Type", size="Significance")
  } else if (v_c == "v") {
    p <- sccomp_df %>%
      ggplot2::ggplot(aes(x = v_effect, y = y, col = diffexpressed, size=v_FDR)) +
      ggplot2::geom_vline(xintercept = 0, col = "gray") +
      ggplot2::geom_point() +
      ggplot2::geom_errorbar(
        aes(xmin = v_lower, xmax = v_upper),
        width = 0.5, position = pd, orientation = "y", size=0.5
      ) +
      ggplot2::scale_y_continuous(
        breaks = sccomp_df |> dplyr::pull(y),
        labels = sccomp_df |> dplyr::pull("cluster_annot")) +
      ggplot2::scale_size_continuous(trans = ggforce::trans_reverser("log10")) +
      ggplot2::xlab("Credible Interval") +
      ggplot2::theme_classic()  +
      ggplot2::scale_color_manual(values = c("dn" = dn, "ns" = "gray", "up" = up)) +
      ggplot2::theme(
        axis.title.x = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        axis.title.y = element_blank() ,
        axis.text.y = element_text(size = 14),
      ) +
      ggplot2::labs(color = "Type", size="Significance")
  }
  
  return (p)
}

clean_hallmark_names <- function(geneset_names) {
  geneset_names <- str_remove_all(geneset_names, pattern = "HALLMARK_")
  geneset_names <- str_replace_all(geneset_names, pattern = "_", replacement = " ")
  geneset_names <- str_to_sentence(geneset_names)
  geneset_names <- str_replace(geneset_names, "Pi3k" , "PI3K") %>%
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
  return(geneset_names)
}

#' Extracts OOB accuracy and Gini-based important features given a dataframe and target column.
#' @param data: A dataframe with feature columns and one column that is the target 
#' @param target_col: A string specifying the name of the target column
#' @param target_class: A string specifying the name of the target class
#'
#' @export
rf_data <- function(data, target_col, target_class) {
  
  # Separate features and target
  features <- data %>% dplyr::select(-target_col)
  data$target <- ifelse(data[[target_col]] == target_class, 1, 0)
  target <- data$target
  
  # Train Random Forest
  rf <- randomForest(x = features, y = as.factor(target), importance = TRUE)
  accuracy <- 1 - rf$err.rate[nrow(rf$err.rate), "OOB"]
  
  # Get feature importances
  importance <- importance(rf)
  
  # Return sorted importances
  importance_df <- data.frame(
    feature = rownames(importance),
    importance = importance[, "MeanDecreaseGini"]
  ) %>%
    arrange(desc(importance))
  
  return(list(acc = accuracy, imps = importance_df))
}
