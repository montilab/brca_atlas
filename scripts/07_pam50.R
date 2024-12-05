library(genefu)
library(ComplexHeatmap)

PATH <- file.path(Sys.getenv("MLAB"), "projects/brcameta/brca_atlas/")
brca_pb <- readRDS(file.path(PATH, "data/pseudobulk/combined_data_avg.rds"))

# Identifying Subtypes
data(pam50)
data(pam50.robust)
pam50_genes <- pam50$centroids.map$probe
pam50_seurat <- brca_pb[pam50_genes,]
pam50_mat <- brca_pb@assays$RNA$data
annot <- pam50$centroids.map
annot <- annot %>% dplyr::rename(Gene.Symbol = probe) %>% dplyr::select(Gene.Symbol, EntrezGene.ID)

pam50_subtypes <- molecular.subtyping(
  sbt.model = "pam50",
  data = t(pam50_mat),  # Transpose the matrix
  annot = annot,
  do.mapping = FALSE
)

pam50_df <- data.frame(donor = names(pam50_subtypes$subtype), pam50 = unname(pam50_subtypes$subtype))
pam50_df$donor <- str_replace_all(pam50_df$donor, "-", "_")
write.csv(pam50_df, file.path(PATH, "data/metadata/pam50.csv"))

# Comparing with author's annotations
all_metadata <- read.csv(file.path(PATH, "data/metadata/all_metadata.csv"))
all_metadata <- all_metadata %>% column_to_rownames("X")
all_metadata$donor <- str_replace_all(all_metadata$donor, "-", "_") 
all_metadata <- all_metadata %>% dplyr::left_join(., y = pam50_df, by = "donor")
write.csv(all_metadata, file.path(PATH, "data/metadata/all_metadata.csv"), row.names = TRUE)
pam50_comparison <- pam50_df %>% dplyr::inner_join(., y = all_metadata, by = "donor", multiple = "first")

cont_mat <- pam50_comparison %>% 
  dplyr::select(subtype_new, pam50.x) %>% 
  table %>%
  as.matrix

col_fun <- colorRamp2(c(min(cont_mat), mean(cont_mat), max(cont_mat)), 
                      c("blue", "white", "red"))
# Create heatmap
heatmap <- Heatmap(cont_mat,
                   name = "Count",
                   row_title = "Author Subtype",
                   column_title = "PAM50 Subtype",
                   col = col_fun,
                   row_order = c("TNBC", "ER+", "HER2+", "Unassigned"),
                   cluster_columns = TRUE,
                   show_row_names = TRUE,
                   show_row_dend = FALSE,
                   show_column_names = TRUE,
                   row_names_side = "left",
                   column_names_side = "top",
                   cell_fun = function(j, i, x, y, width, height, fill) {
                     grid.text(sprintf("%.0f", cont_mat[i, j]), x, y, gp = gpar(fontsize = 10))
                   })

# Draw the heatmap
png(file.path(PATH, "results/annotation/pam50.png"), width=1000, height=1000, res = 200)
draw(heatmap)
dev.off()
