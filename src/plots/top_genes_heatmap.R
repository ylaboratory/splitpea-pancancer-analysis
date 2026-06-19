library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(readr)

cosmic_data <- read_tsv("/ref/Cosmic_MutantCensus_v99_GRCh37.tsv.gz")
unique_driver_genes <- unique(cosmic_data$GENE_SYMBOL)
top_genes_df <- read_csv('blue_heatmap_ordered.csv')
top_genes_df <- as.data.frame(top_genes_df)
gene_names <- top_genes_df$symbol
rownames(top_genes_df) <- gene_names
degree_cols <- grep('norm_degreed', colnames(top_genes_df), value = TRUE)
prop_cols <- gsub('norm_degreed', 'norm_degree_prop', degree_cols)
# # v2
#degree_cols <- grep('presence_prop', colnames(top_genes_df), value = TRUE)
#prop_cols <- gsub('presence_prop', 'presence_prop', degree_cols)
#heatmap_data <- top_genes_df %>% select(all_of(prop_cols))
heatmap_data = top_genes_df[, prop_cols]
colnames(heatmap_data) <- degree_cols
x_labels <- gsub('_norm_degreed', '', degree_cols)
x_labels <- gsub('_', '+', x_labels)
x_labels <- gsub('\\+.*', '', x_labels)   # drop "+" and everything after it
# # v2
#x_labels <- gsub('_presence_prop', '', degree_cols)
#x_labels <- gsub('_', '+', x_labels)
#x_labels <- gsub('\\+.*', '', x_labels)
gene_names <- gsub('\\+.*', '', gene_names)  # same cleanup on row labels
chaos_fraction <- top_genes_df$edge_count_chaos /
  (top_genes_df$edge_count_chaos + top_genes_df$edge_count_neg + top_genes_df$edge_count_pos)
neg_fraction <- top_genes_df$edge_count_neg /
  (top_genes_df$edge_count_chaos + top_genes_df$edge_count_neg + top_genes_df$edge_count_pos)
pos_fraction <- top_genes_df$edge_count_pos /
  (top_genes_df$edge_count_chaos + top_genes_df$edge_count_neg + top_genes_df$edge_count_pos)
fractions <- data.frame(
  Chaos = chaos_fraction,
  Negative = neg_fraction,
  Positive = pos_fraction
)
right_annotation <- rowAnnotation(
  `Edge Types` = anno_barplot(as.matrix(fractions), beside = FALSE,
    gp = gpar(fill = c("#ffd400", "firebrick3", "#6ac072")),
    width = unit(2, "cm"))
)
ht <- Heatmap(as.matrix(heatmap_data),
  name = "proportion",
  col = colorRampPalette(c("#d8e1e8", "#3A527E"))(100),
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.rect(x, y, width, height, gp = gpar(fill = fill, col = "black"))
  },
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 12,
    fontface = ifelse(gene_names %in% unique_driver_genes, "bold", "plain")),
  column_names_side = "bottom",
  column_names_gp = gpar(fontsize = 12),
  column_names_rot = 90,
  column_labels = x_labels,
  row_labels = gene_names,
  right_annotation = right_annotation,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_dend = FALSE,
  show_column_dend = FALSE,
  width = unit(ncol(heatmap_data) * 0.5, "cm"),
  height = unit(nrow(heatmap_data) * 0.5, "cm"))
edge_legend <- Legend(labels = c("Chaos", "Negative", "Positive"),
  title = "Edge Types",
  legend_gp = gpar(fill = c("#ffd400", "firebrick3", "#6ac072")))
draw(ht, annotation_legend_list = list(edge_legend), annotation_legend_side = "right", heatmap_legend_side = "right")