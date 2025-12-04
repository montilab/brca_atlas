library(stringr)

#' Shorten SingleR annotations to broader categories.
#'
#' This function takes a vector of SingleR cell type annotations and simplifies them
#' into more general categories like "Epithelial", "CD4 T-Cells", "Stromal", etc.
#'
#' @param singler_annot A character vector of SingleR cell type annotations.
#' @return A character vector with shortened SingleR annotations.
#' @import stringr
#' @export
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

#' Broadly categorize SingleR annotations into Epithelial, Immune, or Stromal.
#'
#' This function simplifies SingleR cell type annotations into three broad categories:
#' Epithelial, Immune, and Stromal.
#'
#' @param singler_annot A character vector of SingleR cell type annotations.
#' @return A character vector with broad SingleR annotations.
#' @import stringr
#' @export
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

#' Shorten CellTypist annotations to broader categories.
#'
#' This function takes a vector of CellTypist cell type annotations and simplifies them
#' into more general categories like "Epithelial", "CD4 T-Cells", "Stromal", etc.
#'
#' @param celltypist_annot A character vector of CellTypist cell type annotations.
#' @return A character vector with shortened CellTypist annotations.
#' @import stringr
#' @export
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

#' Broadly categorize CellTypist annotations into Epithelial, Immune, or Stromal.
#'
#' This function simplifies CellTypist cell type annotations into three broad categories:
#' Epithelial, Immune, and Stromal.
#'
#' @param celltypist_annot A character vector of CellTypist cell type annotations.
#' @return A character vector with broad CellTypist annotations.
#' @import stringr
#' @export
celltypist_broad <- function(celltypist_annot) {
  # Epi
  celltypist_annot[str_detect(celltypist_annot, "LummHR|basal|Lumsec")] <- "Epithelial"
  # Imm
  celltypist_annot[str_detect(celltypist_annot, "CD4|CD8|DC|GD|NK|bmem|b_naive|plasma|Macro|Mono|mye-prol|T_prol|Mast|Neutrophil")] <- "Immune"
  # Strom
  celltypist_annot[str_detect(celltypist_annot, "Fibro|pericytes|Vas|vsmc|Lymph")] <- "Stromal"
  return(celltypist_annot)
}

#' Shorten Author annotations to broader categories.
#'
#' This function takes a vector of Author cell type annotations and simplifies them
#' into more general categories like "Epithelial", "CD4 T-Cells", "Stromal", etc.
#'
#' @param author_annot A character vector of Author cell type annotations.
#' @return A character vector with shortened Author annotations.
#' @import stringr
#' @export
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

#' Standardize Author annotations.
#'
#' This function takes a vector of Author cell type annotations and standardizes them
#' to consistent terms. Handles NAs.
#'
#' @param author_annot A character vector of Author cell type annotations.
#' @return A character vector with standardized Author annotations.
#' @export
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

#' Broadly categorize Author annotations into Epithelial, Immune, or Stromal.
#'
#' This function simplifies Author cell type annotations into three broad categories:
#' Epithelial, Immune, and Stromal.
#'
#' @param author_annot A character vector of Author cell type annotations.
#' @return A character vector with broad Author annotations.
#' @import stringr
#' @export
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

#' Scale values in a matrix to a specified range.
#'
#' This function scales the values in a matrix to a new range, defaulting to -1 to 1.
#'
#' @param x A numeric matrix.
#' @param new_min The minimum value of the new range (default: -1).
#' @param new_max The maximum value of the new range (default: 1).
#' @return A matrix with values scaled to the new range.
#' @export
scale_to_range <- function(x, new_min = -1, new_max = 1) {
  x_min <- apply(x, 2, min)
  x_max <- apply(x, 2, max)
  t((t(x) - x_min) / (x_max - x_min) * (new_max - new_min) + new_min)
}

#' Prepare data for sccomp plotting.
#'
#' This function prepares a data frame of sccomp results for plotting, filtering by
#' parameter, identifying differentially expressed clusters, and ranking by effect size.
#'
#' @param sccomp_df A data frame containing sccomp results. Must include columns
#'   "parameter", "c_FDR", "c_effect", and an annotation column.
#' @param param The parameter to filter the data frame by.
#' @param v_c Either "c" for cell counts or "v" for variance. Determines which
#'   columns to use for filtering and ranking. Default is "c".
#' @param fdr The false discovery rate threshold for determining differential expression.
#'   Default is 0.05.
#' @param order An optional vector specifying the order of cluster annotations. If provided,
#'   the cluster_annot column will be converted to a factor with the specified levels,
#'   and the data frame will be ordered accordingly.
#'
#' @return A data frame prepared for plotting, with columns for differential expression
#'   status ("diffexpressed") and y-axis position ("y").
#'
#' @importFrom dplyr %>% filter mutate case_when pull
#' @export
sccomp_prep_data <- function(sccomp_df, param, annot_col = "cluster_annot", v_c = "c", fdr = 0.05, order = NULL) {

  stopifnot(all(c("parameter", "c_FDR", "c_effect") %in% colnames(sccomp_df)))

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
    sccomp_df[[annot_col]] <- factor(sccomp_df[[annot_col]], levels = order)
    sccomp_df <- sccomp_df[order(sccomp_df[[annot_col]]),]
  }
  return(sccomp_df)
}

#' Create sccomp plot.
#'
#' This function generates a ggplot object visualizing sccomp results, including
#' effect sizes, credible intervals, and differential expression status.
#'
#' @param sccomp_df A data frame containing sccomp results, prepared using
#'   `sccomp_prep_data`. Must include columns "c_FDR", "c_effect", "diffexpressed",
#'   "c_upper", "c_lower", and "cluster_annot".
#' @param cluster_annot Character, name of annotation column.
#' @param v_c Either "c" for cell counts or "v" for variance. Determines which
#'   columns to use for plotting. Default is "c".
#' @param up The color for upregulated clusters (default: "#FF9999").
#' @param dn The color for downregulated clusters (default: "#56B4E9").
#' @param sort_y A logical value indicating whether to sort the y-axis by effect size
#'   (default: TRUE). If FALSE, clusters will be plotted in the order they appear
#'   in the data frame.
#'
#' @return A ggplot object visualizing the sccomp results.
#'
#' @import ggplot2
#' @import ggforce
#' @importFrom dplyr %>% pull
#' @export
sccomp_plot_data <- function(sccomp_df,
                             annot_col = "cluster_annot",
                             v_c = "c",
                             up = "#FF9999",
                             dn = "#56B4E9",
                             sort_y = TRUE) {

  stopifnot(all(c("c_FDR", "c_effect", "diffexpressed", "c_upper", "c_lower") %in% colnames(sccomp_df)))

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
        labels = sccomp_df |> dplyr::pull(annot_col)) +
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
        labels = sccomp_df |> dplyr::pull(annot_col)) +
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

#' Clean Hallmark gene set names.
#'
#' This function cleans up the names of Hallmark gene sets by removing prefixes,
#' replacing underscores, and converting to sentence case.  It also standardizes
#' some common abbreviations.
#'
#' @param geneset_names A character vector of Hallmark gene set names.
#' @return A character vector of cleaned gene set names.
#' @import stringr
#' @export
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

#' Calculate Cell Type Diversity Statistic
#'
#' This function calculates the cell type diversity statistic for each sample
#' from a matrix with samples specified in rows and cell types specified in
#' the columns. The matrix can be given as an input or generated from a Seurat object
#' or SingleCellExperiment object. If the column sums of the matrix do not equal to 1,
#' the function will return a warning message to fix the matrix.
#'
#'
#' @param normprop A matrix of the normalized cell type proportions with samples by row and
#' cell types by column, or SingleCellExperiment Object or Seurat Object to pull data from.
#' @param sample Variable name for sample ID if extracting from SingleCellExperiment or Seurat Object.
#' @param cell.type Variable name for cell type if extracting variable from SingleCellExperiment or Seurat Object.
#' @param metadata A matrix of the sample level metadata information including sample by row
#' and sample ID variable and other metadata variables by column.
#' @return A matrix of the cell type diversity statistics for each sample joined with metadata
#' information if provided
#' @export
CTDS.score <- function(dataobj,
                       sample = "sample.ID",
                       cell.type = "ct.consensus",
                       metadata = NULL){
  require(tidyverse)
  if(class(dataobj) == "matrix" | class(dataobj) == "table" | class(dataobj) == "data.frame"){
    normprop.table <- dataobj
    #return(normprop.table)
  }else if(class(dataobj) == "SingleCellExperiment"){
    normprop.table <- prop.table(table(dataobj@colData[,sample], dataobj@colData[,cell.type]), margin = 1)
    #return(normprop.table)
  }else if(class(dataobj) == "Seurat"){
    normprop.table <- prop.table(table(dataobj@meta.data[,sample], dataobj@meta.data[,cell.type]), margin = 1)
    #return(normprop.table)
  }else{
    warning("The data object is not a matrix, Seurat object, or a SingleCellExperiment object")
  }

  #check if column sums equal to 1
  normprop.sum <- apply(normprop.table, 1, sum)

  if(all(normprop.sum) == 1){
    message("The normalized proportions for each sample add up to 1")
  }else{
    warning("The normalized proportions for each sample do not add up to 1")
  }

  div.res <- apply(normprop.table, 1, function(x){(-sum(x*log(x), na.rm = T)/log(ncol(normprop.table))-1)})
  if(is.null(metadata) == TRUE){
    div.mat <- tibble::enframe(div.res, name = sample, value = "statistic")
    return(div.mat)
  } else{
    div.meta <- tibble::enframe(div.res, name = sample, value = "statistic") %>%
      dplyr::full_join(metadata, by = sample)
    return(div.meta)
  }


}
#' Adjust P-values for Multiple Linear Models
#'
#' This function takes a list of linear model fits and adjusts the p-values
#' for all coefficients across all models using a specified method.
#'
#' @param model_list A named list of linear model objects (e.g., from lm() or glm()).
#' @param method A character string specifying the adjustment method. Default is "fdr".
#'   Possible values include "holm", "hochberg", "hommel", "bonferroni", "BH", "BY",
#'   "fdr", "none". See ?p.adjust for more details.
#'
#' @return A list of data frames, each containing the coefficient summary for a model,
#'   including the original statistics and an additional column with adjusted p-values.
#'
#' @details This function extracts p-values from all models, adjusts them using
#'   the specified method, and then adds the adjusted p-values back to each
#'   model's coefficient summary. The function preserves the names of the models
#'   and coefficients in the output.
#'
#' @export
p_adjust_list <- function(model_list, method = "fdr") {

  stopifnot(is.list(model_list))

  n_covar <- length(model_list[[1]][["coefficients"]])
  n_models <- length(model_list)

  p_matrix <- matrix(data = NA ,ncol = n_models, nrow = n_covar)
  colnames(p_matrix) <- names(model_list)
  rownames(p_matrix) <- names(model_list[[1]][["coefficients"]])

  for (i in seq_along(model_list)) {
    model_fit <- model_list[[i]]
    coef_df <- summary(model_fit)$coef
    p_vals <- coef_df[,ncol(coef_df)]
    p_matrix[,i] <- p_vals
  }

  p_adj_matrix <- t(apply(p_matrix, 1, p.adjust, method = method, simplify = TRUE))

  coef_dfs <- list()
  for (i in seq_along(model_list)) {
    model_fit <- model_list[[i]]
    model_name <- names(model_list)[[i]]
    coef_df <- summary(model_fit)$coef
    coef_df <- cbind(coef_df, "p.adj" = p_adj_matrix[,i])
    coef_dfs[[model_name]] <- coef_df
  }

  return(coef_dfs)
}

#' Filter Coefficient Data Frame Based on P-value and Direction
#'
#' This function filters a coefficient data frame based on a p-value threshold
#' and the direction of the effect (positive or negative).
#'
#' @param coef_df A data frame containing coefficient information, typically
#'   from a model summary.
#' @param p_column_name Character string specifying the name of the p-value column.
#'   Default is "p.adj".
#' @param p_threshold Numeric value specifying the p-value threshold for filtering.
#'   Default is 0.05.
#' @param pos_only Logical indicating whether to keep only positive effects.
#'   Default is TRUE.
#' @param neg_only Logical indicating whether to keep only negative effects.
#'   Default is FALSE.
#'
#' @return The filtered coefficient data frame if conditions are met, otherwise NULL.
#'
#' @details
#' The function first checks if the last row's p-value is above the threshold.
#' If so, it returns NULL. Then, based on the pos_only and neg_only parameters,
#' it filters for positive or negative effects. If both pos_only and neg_only
#' are FALSE, it returns the data frame without directional filtering.
#'
#' @examples
#' coef_df <- data.frame(
#'   Estimate = c(0.5, -0.3, 0.7),
#'   `Pr(>|t|)` = c(0.01, 0.06, 0.03),
#'   p.adj = c(0.02, 0.08, 0.04)
#' )
#'
#' # Filter for positive effects with adjusted p-value < 0.05
#' filter_fit(coef_df, p_column_name = "p.adj", p_threshold = 0.05, pos_only = TRUE)
#'
#' # Filter for negative effects with p-value < 0.1
#' filter_fit(coef_df, p_column_name = "Pr(>|t|)", p_threshold = 0.1, pos_only = FALSE, neg_only = TRUE)
#'
#' @export
filter_fit <- function(coef_df,
                       p_column_name = "p.adj",
                       p_threshold = 0.05,
                       pos_only = TRUE,
                       neg_only = FALSE) {

  if (coef_df[nrow(coef_df),p_column_name] > p_threshold) {
    return(NULL)
  }

  if(pos_only) {
    if (coef_df[nrow(coef_df),1] < 0) {
      return(NULL)
    } else {
      return(coef_df)
    }
  }
  if(neg_only) {
    if (coef_df[nrow(coef_df),1] > 0) {
      return(NULL)
    } else {
      return(coef_df)
    }
  }
}

fishers_meta_p <- function(pvals) {
  # pvals: a vector of p-values, e.g. c(0.01, 0.03, 0.25)
  X <- -2 * sum(log(pvals))
  df <- 2 * length(pvals)
  return(pchisq(X, df=df, lower.tail = FALSE))
}
