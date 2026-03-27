library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
library(ggridges)

data <- read.csv("/results/edge_stats.csv", stringsAsFactors = TRUE)

filtered_data <- data %>% 
  filter((Positive_Count + Negative_Count) / Total_Samples >= 0.5)

middle_bin_data <- filtered_data %>%
  filter(Stability > -0.25, Stability <= 0.25) %>%
  mutate(gene_list = str_extract_all(Edge_Symbol, "[A-Z0-9]+")) %>%
  unnest(gene_list) %>%
  rename(gene = gene_list)


gene_cancer_counts <- middle_bin_data %>%
  distinct(gene, cancer) %>%
  count(gene, name = "num_cancers")

genes_unique_to_one_cancer <- gene_cancer_counts %>%
  filter(num_cancers == 1) %>%
  select(gene)

unique_gene_counts <- unique_gene_counts %>%
  arrange(desc(unique_gene_count)) %>%
  mutate(cancer = factor(cancer, levels = rev(cancer))) 

unique_gene_counts <- unique_gene_counts %>%
  arrange(desc(unique_gene_count)) %>%
  mutate(cancer = factor(cancer, levels = rev(cancer)))
  
r_grey <- "#7C7E7E"

p_unique_middle_bin <- ggplot(
  unique_gene_counts,
  aes(x = cancer, y = unique_gene_count)
) +
  geom_col(fill = r_grey) +                
  labs(
    y = "Number of Genes From Non-Polarized Edges Unique to This Cancer",
    x = NULL
  ) +
  coord_flip() +
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10)
  )

p_unique_middle_bin


nonunique_genes <- gene_cancer_counts %>%
  filter(num_cancers > 1)

ridge_df <- middle_bin_data %>%
  distinct(gene, cancer) %>%
  inner_join(nonunique_genes, by = "gene") %>%
  mutate(
    n_other_cancers = num_cancers - 1,
    cancer = factor(
      cancer,
      levels = levels(unique_gene_counts$cancer)
    )
  )

p_ridge_nonunique <- ggplot(
  ridge_df,
  aes(x = n_other_cancers, y = cancer)
) +
  geom_density_ridges(
    stat = "binline",
    scale = 1,
    bins = 10,
    fill = r_grey,   
    color = "black",  
  ) +
  labs(
    x = "Number of other cancers each non-unique gene appears in",
    y = NULL
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    axis.text.y = element_text(size = 10)
  )

p_ridge_nonunique
