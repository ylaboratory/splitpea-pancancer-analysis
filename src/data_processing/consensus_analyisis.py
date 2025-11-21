import os 
import pandas as pd
import pickle 

def get_folders(directory):
    folders_list = []

    folders = [f for f in os.listdir(directory) if os.path.isdir(os.path.join(directory, f))]

    for folder in folders:
        subfolders = [f for f in os.listdir(os.path.join(directory, folder)) if os.path.isdir(os.path.join(directory, folder, f))]
        for subfolder in subfolders:
            folders_list.append((folder, subfolder))

    return folders_list


directory_path = "/results/splitpea_outputs/"
folders_tuples = get_folders(directory_path)
print(folders_tuples)

file_path = '/ref/hsa_mapping_all.txt'
gene_mappings = pd.read_csv(file_path, sep = '\t')
symbol_to_entrez = dict(zip(gene_mappings['symbol'], gene_mappings['entrez']))
entrez_to_symbol = {value: key for key, value in symbol_to_entrez.items() if not pd.isna(value)}

cancers = folders_tuples[0]

base_dir = "/results/splitpea_outputs/" 

for cancers in folders_tuples:
    
    cancer_folder = cancers[0]
    subfolder = cancers[1]
    
    neg_path = os.path.join(base_dir, cancer_folder, subfolder, "consensus_network_neg.pickle")
    pos_path = os.path.join(base_dir, cancer_folder, subfolder, "consensus_network_pos.pickle")
    
    print(f"Loading {cancer_folder} consensus networks...")
    
    with open(neg_path, 'rb') as f_neg, open(pos_path, 'rb') as f_pos:
        print("opening...")
        c_consensus_neg = pickle.load(f_neg)
        c_consensus_pos = pickle.load(f_pos)

edge_counts = {}
cancer_num_graphs = {}  

for cancers in folders_tuples:
    cancer_folder = cancers[0]
    subfolder = cancers[1]

    neg_path = os.path.join(base_dir, cancer_folder, subfolder, "consensus_network_neg.pickle")
    pos_path = os.path.join(base_dir, cancer_folder, subfolder, "consensus_network_pos.pickle")

    print(f"Loading {cancer_folder} consensus networks...")

    with open(neg_path, 'rb') as f_neg, open(pos_path, 'rb') as f_pos:
        print("Opening consensus networks...")
        c_consensus_neg = pickle.load(f_neg)
        c_consensus_pos = pickle.load(f_pos)
    
    num_graphs_neg = c_consensus_neg.graph.get('num_graphs', 1)
    num_graphs_pos = c_consensus_pos.graph.get('num_graphs', 1)
    
    if num_graphs_neg != num_graphs_pos:
        print(f"Warning: 'num_graphs' differs between positive and negative networks for {cancer_folder}")
    num_graphs = max(num_graphs_neg, num_graphs_pos)
    cancer_num_graphs[cancer_folder] = num_graphs  

    print("Processing negative edges...")
    for edge in c_consensus_neg.edges():
        num_neg = c_consensus_neg.edges[edge].get('num_neg', 0)
        if edge not in edge_counts:
            edge_counts[edge] = {}
        if cancer_folder not in edge_counts[edge]:
            edge_counts[edge][cancer_folder] = {'Positive_Count': 0, 'Negative_Count': 0}
        edge_counts[edge][cancer_folder]['Negative_Count'] += num_neg
    
    print("Processing positive edges...")
    for edge in c_consensus_pos.edges():
        num_pos = c_consensus_pos.edges[edge].get('num_pos', 0)
        if edge not in edge_counts:
            edge_counts[edge] = {}
        if cancer_folder not in edge_counts[edge]:
            edge_counts[edge][cancer_folder] = {'Positive_Count': 0, 'Negative_Count': 0}
        edge_counts[edge][cancer_folder]['Positive_Count'] += num_pos

data = []
for edge, cancers in edge_counts.items():
    for cancer, counts in cancers.items():
        data.append({
            'Edge': edge,
            'Cancer': cancer,
            'Positive_Count': counts['Positive_Count'],
            'Negative_Count': counts['Negative_Count'],
            'Total_Samples': cancer_num_graphs.get(cancer, 1) 
        })

df = pd.DataFrame(data)

df['Raw_Difference'] = df['Positive_Count'] - df['Negative_Count']
df['Normalized_Difference'] = df['Raw_Difference'] / df['Total_Samples']
df['Rank Pos'] = df.groupby('Cancer')['Normalized_Difference'].rank(method='first', ascending=False)
df['Rank Neg'] = df.groupby('Cancer')['Normalized_Difference'].rank(method='first', ascending=True)

df_analysis = df.copy()

edge_totals = df_analysis.groupby('Edge').agg({'Positive_Count': 'sum', 'Negative_Count': 'sum'}).reset_index()
edge_totals.rename(columns={'Positive_Count': 'Total_Positive_Count', 'Negative_Count': 'Total_Negative_Count'}, inplace=True)

df_analysis = df_analysis.merge(edge_totals, on='Edge')

df_analysis['Positive_Specificity_Score'] = df_analysis['Positive_Count'] / df_analysis['Total_Positive_Count']
df_analysis['Negative_Specificity_Score'] = df_analysis['Negative_Count'] / df_analysis['Total_Negative_Count']

df_analysis['Positive_Specificity_Score'] = df_analysis['Positive_Specificity_Score'].fillna(0)
df_analysis['Negative_Specificity_Score'] = df_analysis['Negative_Specificity_Score'].fillna(0)

df_analysis['Positive_Rate'] = df_analysis['Positive_Count'] / df_analysis['Total_Samples']
df_analysis['Negative_Rate'] = df_analysis['Negative_Count'] / df_analysis['Total_Samples']

df_analysis['Adjusted_Positive_Specificity_Score'] = df_analysis['Positive_Specificity_Score'] * df_analysis['Positive_Rate']
df_analysis['Adjusted_Negative_Specificity_Score'] = df_analysis['Negative_Specificity_Score'] * df_analysis['Negative_Rate']

df_analysis['Net_Edge_Score'] = (df_analysis['Positive_Count'] - df_analysis['Negative_Count']) / df_analysis['Total_Samples']
df_analysis['Net_Specificity_Score'] = df_analysis['Positive_Specificity_Score'] - df_analysis['Negative_Specificity_Score']
df_analysis['Adjusted_Net_Specificity_Score'] = df_analysis['Adjusted_Positive_Specificity_Score'] - df_analysis['Adjusted_Negative_Specificity_Score']

df_analysis['Rank Neg Spec'] = df_analysis['Adjusted_Net_Specificity_Score'].rank(method='first', ascending=True)
df_analysis['Rank Pos Spec'] = df_analysis['Adjusted_Net_Specificity_Score'].rank(method='first', ascending=False)

df_analysis["Stability"] = df_analysis["Positive_Rate"] - df_analysis["Negative_Rate"]

df_analysis.to_csv("/results/edge_stats.csv")

"""
Edge
Cancer
Positive_Count: # of samples with that edge in that cancer 
Negative_Count
Total_Samples: # of samples in the cancer 

Raw_Difference:  (Positive_Count - Negative_Count).
Normalized_Difference: (Raw_Difference / Total_Samples).
Total_Positive_Count: # samples edge was positive across all cancers.
Total_Negative_Count: same with neg

Positive_Specificity_Score: (Positive_Count / Total_Positive_Count).
Negative_Specificity_Score: (Negative_Count / Total_Negative_Count).
Positive_Rate: (Positive_Count / Total_Samples).
Negative_Rate:  (Negative_Count / Total_Samples).

Adjusted_Positive_Specificity_Score:  (Positive_Specificity_Score * Positive_Rate).
Adjusted_Negative_Specificity_Score:  (Negative_Specificity_Score * Negative_Rate).
Net_Edge_Score:  ((Positive_Count - Negative_Count) / Total_Samples).
Net_Specificity_Score: (Positive_Specificity_Score - Negative_Specificity_Score).
Adjusted_Net_Specificity_Score:  (Adjusted_Positive_Specificity_Score - Adjusted_Negative_Specificity_Score).

Stability: Postive_Rate - Negative_Rate

"""


