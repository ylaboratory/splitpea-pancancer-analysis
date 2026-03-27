library(dplyr)
library(stringr)
library(tidyr)
library(clusterProfiler)
library(org.Hs.eg.db)

extract_genes <- function(edge_str) {
  unique(unlist(strsplit(gsub("[()']", "", edge_str), ",\\s*")))
}

run_enrichGO <- function(genes, background) {
  enrichGO(
    gene          = genes,
    universe      = background,
    OrgDb         = org.Hs.eg.db,
    keyType       = "ENTREZID",
    ont           = "BP",
    pAdjustMethod = "BH",
    qvalueCutoff  = 0.05,
    readable      = TRUE
  )
}

counts_per_edge_neg <- read.csv('/data/counts_per_edge_neg.csv',
                                stringsAsFactors = FALSE) # from /plots/edge_conservation.py

bg_genes_neg <- unique(unlist(lapply(counts_per_edge_neg$Edge, extract_genes)))
cat("Neg background - unique genes:", length(bg_genes_neg), "\n")

case1_neg <- counts_per_edge_neg %>% filter(Threshold == 0, Cancer_Count == 28)
genes_case1_neg <- unique(unlist(lapply(case1_neg$Edge, extract_genes)))
cat("Case 1 neg - edges:", nrow(case1_neg), "| genes:", length(genes_case1_neg), "\n")

case2_neg <- counts_per_edge_neg %>% filter(Cancer_Count >= 10, Threshold >= 0.75)
genes_case2_neg <- unique(unlist(lapply(case2_neg$Edge, extract_genes)))
cat("Case 2 neg - edges:", nrow(case2_neg), "| genes:", length(genes_case2_neg), "\n")

ego_case1_neg <- run_enrichGO(genes_case1_neg, bg_genes_neg)
ego_case2_neg <- run_enrichGO(genes_case2_neg, bg_genes_neg)

counts_per_edge_pos <- read.csv('/data/counts_per_edge_pos.csv',
                                stringsAsFactors = FALSE)

bg_genes_pos <- unique(unlist(lapply(counts_per_edge_pos$Edge, extract_genes)))
cat("Pos background - unique genes:", length(bg_genes_pos), "\n")

case1_pos <- counts_per_edge_pos %>% filter(Threshold == 0, Cancer_Count == 28)
genes_case1_pos <- unique(unlist(lapply(case1_pos$Edge, extract_genes)))
cat("Case 1 pos - edges:", nrow(case1_pos), "| genes:", length(genes_case1_pos), "\n")

case2_pos <- counts_per_edge_pos %>% filter(Cancer_Count >= 10, Threshold >= 0.75)
genes_case2_pos <- unique(unlist(lapply(case2_pos$Edge, extract_genes)))
cat("Case 2 pos - edges:", nrow(case2_pos), "| genes:", length(genes_case2_pos), "\n")

ego_case1_pos <- run_enrichGO(genes_case1_pos, bg_genes_pos)
ego_case2_pos <- run_enrichGO(genes_case2_pos, bg_genes_pos)

if (!dir.exists("data/go")) dir.create("data/go", recursive = TRUE)

write.csv(ego_case1_neg@result, "data/go/thresh0_all_cancers_loss.csv", row.names = FALSE)
write.csv(ego_case2_neg@result, "data/go/thresh75_10_loss_cancers.csv", row.names = FALSE)
write.csv(ego_case2_pos@result, "data/go/thresh75_10_gain_cancers.csv", row.names = FALSE)

edge_stats <- read.csv("/data/edge_stats.csv",
                       stringsAsFactors = TRUE)

# stability is the same as consistency 
middle_bin_data <- edge_stats %>%
  filter((Positive_Count + Negative_Count) / Total_Samples >= 0.5,
         Stability > -0.25, Stability <= 0.25) %>%
  mutate(gene_list = str_extract_all(Edge_Symbol, "[A-Z0-9]+")) %>%
  unnest(gene_list) %>%
  rename(gene = gene_list) %>%
  dplyr::select(cancer, gene) %>%
  distinct()

bg_entrez <- unique(unlist(lapply(edge_stats$Edge, extract_genes)))

symbol2entrez <- bitr(
  unique(middle_bin_data$gene),
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = org.Hs.eg.db
)

map_to_entrez <- function(symbols) {
  unique(symbol2entrez$ENTREZID[symbol2entrez$SYMBOL %in% symbols])
}

out_dir <- "data/go/cancer_nonpolar"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

for (ca in sort(unique(middle_bin_data$cancer))) {
  cat("Running GO for:", ca, "\n")

  genes_entrez <- map_to_entrez(middle_bin_data$gene[middle_bin_data$cancer == ca])

  if (length(genes_entrez) == 0) {
    cat("  No ENTREZ IDs – skipping.\n")
    next
  }

  ego <- run_enrichGO(genes_entrez, bg_entrez)
  ego_df <- as.data.frame(ego)

  if (nrow(ego_df) == 0) {
    cat("  No significant GO terms – skipping.\n")
    next
  }

  out_file <- file.path(out_dir, paste0("GO_BP_", str_replace_all(ca, "[^A-Za-z0-9]+", "_"), ".csv"))
  write.csv(ego_df, out_file, row.names = FALSE)
  cat("  Saved:", out_file, "\n")
}

consensus_rare_neg <- counts_per_edge_neg %>% filter(Threshold == 1, Cancer_Count <= 2)
consensus_rare_pos <- counts_per_edge_pos %>% filter(Threshold == 1, Cancer_Count <= 2)

genes_consensus_rare <- unique(c(
  unlist(lapply(consensus_rare_neg$Edge, extract_genes)),
  unlist(lapply(consensus_rare_pos$Edge, extract_genes))
))
cat("Consensus rare (thresh=1, <=2 cancers) - genes:", length(genes_consensus_rare), "\n")

ego_consensus_rare <- run_enrichGO(genes_consensus_rare, bg_entrez)
write.csv(as.data.frame(ego_consensus_rare), "data/go/consensus_thresh1_max2_cancers.csv", row.names = FALSE)

polar_data <- edge_stats %>%
  filter((Positive_Count + Negative_Count) / Total_Samples >= 0.5,
         Stability < -0.3 | Stability > 0.3) %>%
  mutate(gene_list = str_extract_all(Edge_Symbol, "[A-Z0-9]+")) %>%
  unnest(gene_list) %>%
  rename(gene = gene_list) %>%
  dplyr::select(cancer, gene) %>%
  distinct()

symbol2entrez_polar <- bitr(
  unique(polar_data$gene),
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = org.Hs.eg.db
)

map_to_entrez_polar <- function(symbols) {
  unique(symbol2entrez_polar$ENTREZID[symbol2entrez_polar$SYMBOL %in% symbols])
}

out_dir_polar <- "data/go/cancer_polar"
if (!dir.exists(out_dir_polar)) dir.create(out_dir_polar, recursive = TRUE)

for (ca in sort(unique(polar_data$cancer))) {
  cat("Running GO for polar edges:", ca, "\n")

  genes_entrez <- map_to_entrez_polar(polar_data$gene[polar_data$cancer == ca])

  if (length(genes_entrez) == 0) {
    cat("  No ENTREZ IDs – skipping.\n")
    next
  }

  ego <- run_enrichGO(genes_entrez, bg_entrez)
  ego_df <- as.data.frame(ego)

  if (nrow(ego_df) == 0) {
    cat("  No significant GO terms – skipping.\n")
    next
  }

  out_file <- file.path(out_dir_polar, paste0("GO_BP_", str_replace_all(ca, "[^A-Za-z0-9]+", "_"), ".csv"))
  write.csv(ego_df, out_file, row.names = FALSE)
  cat("  Saved:", out_file, "\n")
}
