import seaborn as sns
import matplotlib.pyplot as plt
from collections import Counter
import pandas as pd
from scipy.cluster.hierarchy import linkage, leaves_list
from matplotlib.patches import Patch

cell_width = 1   
cell_height = 0.2 
fontsize = 12       

df_analysis = pd.read_csv("/results/edge_stats.csv")  # from consensus_analysis.py

low_thresh = 0.7

pos_rate_df = df_analysis[
    (df_analysis['Positive_Rate'] >= low_thresh)
]

neg_rate_df = df_analysis[
    (df_analysis['Negative_Rate'] >= low_thresh)
]

cancer_thresh = 3 

pos_edge_cancer_count = pos_rate_df.groupby('Edge')['Cancer'].nunique()
neg_edge_cancer_count = neg_rate_df.groupby('Edge')['Cancer'].nunique()

pos_edges_high_cancer = pos_edge_cancer_count[pos_edge_cancer_count >= cancer_thresh]
neg_edges_high_cancer = neg_edge_cancer_count[neg_edge_cancer_count >= cancer_thresh]

common_edges = pos_edges_high_cancer.index.intersection(neg_edges_high_cancer.index)

common_edge_counts = pd.DataFrame({
    'Positive_Cancer_Count': pos_edges_high_cancer.loc[common_edges],
    'Negative_Cancer_Count': neg_edges_high_cancer.loc[common_edges]
})

common_edge_counts = common_edge_counts.reset_index()

gene_counter = Counter()
for edge in common_edge_counts['Edge_Symbol']:
    gene_counter.update(edge)


significant_genes = [gene for gene, count in gene_counter.items() if count >= 10]


color_palette = sns.color_palette("husl", len(significant_genes))
color_map = {}
for gene in gene_counter:
    if gene in significant_genes:
        color_map[gene] = color_palette[significant_genes.index(gene)]
    else:
        color_map[gene] = "grey"


filtered_df_analysis = df_analysis[df_analysis['Edge'].isin(common_edge_counts["Edge"])]
heatmap_data = filtered_df_analysis.pivot(
    index="Edge_Symbol", 
    columns="Cancer", 
    values="Normalized_Difference"
).fillna(0)


row_colors = []
edge_colors = {}
for edge in heatmap_data.index:
    gene1, gene2 = edge
    if gene1 in significant_genes:
        color = color_map[gene1]
    elif gene2 in significant_genes:
        color = color_map[gene2]
    else:
        color = "grey"
    row_colors.append(color)
    edge_colors[edge] = color 


heatmap_data_reset = heatmap_data.reset_index()
heatmap_data_reset['color'] = row_colors

color_counts = {}
unique_colors = set(row_colors)
for color in unique_colors:
    genes_in_color = [g for g, clr in color_map.items() if clr == color]
    total_count = sum(gene_counter[gene] for gene in genes_in_color)
    color_counts[color] = total_count

ordered_colors = [c for c, cnt in sorted(color_counts.items(), key=lambda x: x[1], reverse=True)]
if 'grey' in ordered_colors:
    ordered_colors.remove('grey')
    ordered_colors.append('grey')

print("Ordered Colors:", ordered_colors)  # Debugging


ordered_indices = []
for color in ordered_colors:
    group_df = heatmap_data_reset[heatmap_data_reset['color'] == color]
    if len(group_df) == 0:
        continue
    group_data = group_df.drop(['Edge_Symbol', 'color'], axis=1)
    if len(group_data) > 1:
        Z = linkage(group_data, method='average', metric='euclidean')
        leaves = leaves_list(Z)
        group_order = group_df.index[leaves]
    else:
        group_order = group_df.index
    ordered_indices.extend(group_order)

if isinstance(ordered_indices[0], pd.Index):
    ordered_indices = [idx for sublist in ordered_indices for idx in sublist]

heatmap_data_ordered = heatmap_data.iloc[ordered_indices]
row_colors_ordered = [row_colors[i] for i in ordered_indices]

edges_other = [
    edge for edge in heatmap_data_ordered.index 
    if edge[0] not in significant_genes and edge[1] not in significant_genes
]
edges_sig = [edge for edge in heatmap_data_ordered.index if edge not in edges_other]

ordered_indices_sig = []
ordered_indices_other = []
row_colors_sig = []
row_colors_other = []

for idx, edge in enumerate(heatmap_data_ordered.index):
    color = row_colors_ordered[idx]
    if edge in edges_sig:
        ordered_indices_sig.append(idx)
        row_colors_sig.append(color)
    else:
        ordered_indices_other.append(idx)
        row_colors_other.append(color)

heatmap_data_sig_ordered = heatmap_data_ordered.iloc[ordered_indices_sig]
heatmap_data_other_ordered = heatmap_data_ordered.iloc[ordered_indices_other]

def format_labels(edges, significant_genes):
    new_labels = []
    for edge in edges:
        gene1, gene2 = edge
        if gene1 in significant_genes and gene2 not in significant_genes:
            label = gene2
        elif gene2 in significant_genes and gene1 not in significant_genes:
            label = gene1
        elif gene1 in significant_genes and gene2 in significant_genes:
            label = f"{gene1},{gene2}"
        else:
            label = f"{gene1},{gene2}"
        label = label.replace("'", "").replace("(", "").replace(")", "")
        new_labels.append(label)
    return new_labels

new_labels_sig = format_labels(heatmap_data_sig_ordered.index, significant_genes)
new_labels_other = format_labels(heatmap_data_other_ordered.index, significant_genes)


valid_colors_sig = all(isinstance(color, (tuple, list, str)) for color in row_colors_sig)
valid_colors_other = all(isinstance(color, (tuple, list, str)) for color in row_colors_other)

if not valid_colors_sig:
    print("Error: Invalid colors detected in row_colors_sig.")
if not valid_colors_other:
    print("Error: Invalid colors detected in row_colors_other.")

clustermap_sig = sns.clustermap(
    heatmap_data_sig_ordered,
    cmap="coolwarm",
    row_colors=row_colors_sig,
    yticklabels=new_labels_sig,
    row_cluster=False,
    col_cluster=True,
    tree_kws={"linewidths": 0.},
    figsize=(heatmap_data_sig_ordered.shape[1]*cell_width, 
             heatmap_data_sig_ordered.shape[0]*cell_height)
)

clustermap_sig.ax_heatmap.set_yticklabels(
    clustermap_sig.ax_heatmap.get_ymajorticklabels(),
    fontsize=fontsize
)
clustermap_sig.ax_heatmap.tick_params(axis='y', which='both', right=False)

sig_col_order = clustermap_sig.data2d.columns


heatmap_data_other_ordered = heatmap_data_other_ordered[sig_col_order]

clustermap_other = sns.clustermap(
    heatmap_data_other_ordered,
    cmap="coolwarm",
    row_colors=row_colors_other,
    yticklabels=new_labels_other,
    row_cluster=False,
    col_cluster=False,  
    tree_kws={"linewidths": 0.},
    figsize=(heatmap_data_other_ordered.shape[1]*cell_width, 
             heatmap_data_other_ordered.shape[0]*cell_height)
)

clustermap_other.ax_heatmap.set_yticklabels(
    clustermap_other.ax_heatmap.get_ymajorticklabels(),
    fontsize=fontsize
)
clustermap_other.ax_heatmap.tick_params(axis='y', which='both', right=False)


# Legend for significant genes heatmap
color_gene_map_sig = {}
for gene, color in color_map.items():
    if gene in significant_genes:
        color_gene_map_sig.setdefault(color, []).append(gene)

# Include 'Other Genes' in grey
color_gene_map_sig['grey'] = ['Other Genes']

handles_sig = []
for color in ordered_colors:
    genes = color_gene_map_sig.get(color, [])
    if color != 'grey' and genes:
        label = f"{', '.join(genes)}"
        handles_sig.append(Patch(facecolor=color, edgecolor='black', label=label))
    elif color == 'grey':
        label = 'Other Genes'
        handles_sig.append(Patch(facecolor=color, edgecolor='black', label=label))

clustermap_sig.ax_heatmap.legend(
    handles=handles_sig,
    loc='upper right',
    bbox_to_anchor=(1.5, 1),
    borderaxespad=0,
    fontsize=12,
    title='Gene Legend'
)

# Legend for other genes heatmap
handles_other = [Patch(facecolor="grey", edgecolor='black', label='Other Genes')]
clustermap_other.ax_heatmap.legend(
    handles=handles_other,
    loc='upper right',
    bbox_to_anchor=(1.5, 1),
    borderaxespad=0,
    fontsize=12,
    title='Gene Legend'
)

plt.show()