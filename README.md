# Alternative splicing reshapes protein interaction networks through structured, directional rewiring across cancers

## About
Alternative splicing generates extensive transcript diversity in cancer, but whether the resulting changes to protein-protein interaction (PPI) networks are structured across tumors or idiosyncratic to individual patients remains unresolved. Existing analyses catalog splicing events and their effects on individual interactions, but do not reconstruct patient-specific interaction networks at scale, leaving the systems-level architecture of splicing-driven rewiring uncharacterized. Here, focusing on exon skipping, the most prevalent form of alternative splicing in humans, we analyze 7,949 tumors across 28 cancer types to reconstruct patient-specific rewired PPI networks relative to matched normal tissues, distinguishing interaction gains from losses. We identify widespread and structured remodeling of protein connectivity, including a unique pattern of interaction gains in brain cancers, a conserved pan-cancer axis of interaction loss that converges on GTPase and Ras signaling, and a subset of genes that function as bidirectional switches, where rewiring direction varies across cancer types. At the patient level, tumors exhibiting extreme rewiring, especially pronounced interaction gains, are associated with significantly poorer survival. These results delineate systems-level patterns through which alternative splicing reshapes protein interaction networks and contributes to tumor heterogeneity across cancers.

## Getting started

```bash
conda env create -f env.yml
conda activate splitpea-pancancer
```

## src/run_splitpea

`run_splitpea.sh` runs the full splitpea pipeline for one tissue. Install splitpea from [ylaboratory/splitpea](https://github.com/ylaboratory/splitpea). TCGA and GTEx alternative splicing data can be downloaded from [IRIS](https://github.com/Xinglab/IRIS). COSMIC data can be downloaded from [COSMIC](https://cancer.sanger.ac.uk/cosmic/login) (requires registration).

## src/data_processing

| Script | Description |
|---|---|
| `general_analyisis.py` | Loads sample-level splitpea networks; outputs edge/gene summaries, scatter data, top gene degrees, driver edge info, and heatmap data |
| `consensus_and_consistency_analyisis.py` | Aggregates edge counts across cancers and computes per-edge consistency and other metrics |
| `mutations_mc3.R` | Loads TCGA MC3 mutation data and saves per-sample mutation calls |
| `mutations_drivers.R` | Compares node degrees between mutated drivers, non-mutated drivers, and non-drivers with statistical tests |
| `go_analysis.R` | GO biological process enrichment analysis |

## src/plots

| Script | Description |
|---|---|
| `exon_skipping_dist_plot.R` | Skipped exon event counts per sample across cancer types with per-cohort medians |
| `rewiring_basic_summary_plot.R` | Edge gain/loss/chaos boxplots per cancer type and scatter plots of splicing events with rewiring |
| `stacked_driver_barplots.R` | Proportions of driver-associated edge categories per cancer type |
| `top_genes_heatmap.R` | Heatmap of top genes across cancers with clustering |
| `edge_conservation.py` | Edge conservation across cancers at varying positive/negative rate thresholds |
| `polarizing_edges.py` | Clustermap of edges consistently gained or lost across cancers |
| `nonpolarizing_edges.R` | Genes from non-polarized edges (near-zero consistency) unique to each cancer |
| `plot_outlier_analysis.py` | Analysis and plots of outlier patients with extreme network rewiring|
