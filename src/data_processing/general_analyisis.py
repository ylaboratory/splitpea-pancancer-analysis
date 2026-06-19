import os 
import pandas as pd
import networkx as nx
import pickle
import time
import glob
import sys
from collections import defaultdict
import networkx as nx


def get_files(folder_path, compare):
    file_paths = set()
    file_ind = set()

    for root, dirs, files in os.walk(folder_path):
        for file in files:
            if file not in file_ind and compare in file:
                file_path = os.path.join(root, file)
                file_paths.add(file_path)
                file_ind.add(file)

    return file_paths

directory_path = "/data/splitpea_inputs/"
se_cancer_files = get_files(directory_path, "TCGA")
print(se_cancer_files)
print(len(se_cancer_files))

se_tissue_files = get_files(directory_path, "GTEx")
print(se_tissue_files)
print(len(se_tissue_files))

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

non_nan_counts_dfs = []

for cancer_file in list(se_cancer_files):
    print(cancer_file)
    cancer_name = os.path.basename(cancer_file).split('TCGA_')[1].split('_T')[0]
    df = pd.read_csv(cancer_file, sep='\t')
    
    tcga_columns = df.filter(like='TCGA')
    
    non_nan_counts = tcga_columns.notna().sum()
    
    non_nan_counts_df = pd.DataFrame(non_nan_counts, columns=['non_nan_count'])
    
    non_nan_counts_df['cancer_name'] = cancer_name
    
    non_nan_counts_dfs.append(non_nan_counts_df)


cancer_alt_splice = pd.concat(non_nan_counts_dfs, ignore_index=True)
cancer_alt_splice = cancer_alt_splice.rename(columns={"cancer_name": "Cancer.Types", "non_nan_count": "cancer_tissue_splicemb"})
cancer_alt_splice

cancer_alt_splice.to_csv('/data/cancer_alt_splice.csv', index=False)

row_counts = []
non_nan_counts_dfs = []

for ct in folders_tuples:
    print(ct[0])
    path = "/splitpea_results/"+  '/'+ ct[0] + '/' + ct[1] + '/psis'
    file_list = os.listdir(path)
    print(len(file_list))
    x = 0 
    for file_name in file_list:
        if x% 100 == 0:
            print(x)
        if file_name.endswith('.txt'):
            file_path = os.path.join(path, file_name)
            df = pd.read_csv(file_path, delimiter='\t') 
            filtered_df = df[(df['delta.psi'].abs() > 0.05) & (df['pval'] < 0.05)]
            #row_counts.append((file_name, len(df)))
            row_counts.append((ct[0], len(filtered_df)))
            non_nan_counts_dfs.append((file_name.split('-psi.txt')[0], len(filtered_df)))
        x+=1

cancer_tissue_alt_splice = pd.DataFrame(row_counts, columns=['Cancer.Types', 'cancer_tissue_splice_correctedmb'])
cancer_tissue_alt_splice = cancer_tissue_alt_splice[cancer_tissue_alt_splice.columns[::-1]]

cancer_tissue_alt_splice.to_csv('/data/cancer_tissue_alt_splice.csv', index=False)

edge_df = pd.read_csv('/results/tcga_splitpea_network.largest_ccs.txt',
                      sep="\t")
sample_alt_splice = pd.DataFrame(non_nan_counts_dfs, columns=['sample', 'alternative_se_sig_pd'])

scatter_df = pd.merge(sample_alt_splice, edge_df, on='sample', how='inner')
scatter_df['cancer_type'] = scatter_df['cancer'].str.split('_').str[0]

scatter_df.to_csv('/results/scatter_df.csv', index=False)

base_dir3 = "/ref/hsa_mapping_all.txt"
entrez_map = pd.read_csv(base_dir3, delimiter='\t')

edge_count_neg = defaultdict(int) 
edge_count_pos = defaultdict(int) 
edge_count_chaos = defaultdict(int) 
edge_count_total = defaultdict(int) 
gene_presence_count = defaultdict(int) 
gene_degree_count = defaultdict(int) 

metrics_per_cancer = defaultdict(dict)
    
### fix     
for can in folders_tuples:
    base_dir = "/splitpea/" + can[0] + '/' + can[1] + "/output/"
    base_dir2 = "/splitpea/" + can[0] + '/' + can[1] + "/reference/human_ppi_ddi_bg.pickle"
    

    if len(sys.argv) > 1:
        net_dir = base_dir
        if not os.path.exists(net_dir):
            sys.exit('Error: Input directory containing sample-level networks does not exist.')
    else:
        sys.exit('Error.')

    print('[' + time.strftime("%H:%M:%S", time.localtime()) + "] Loading networks... " + can[0])
    
    c_edge_count_neg = defaultdict(int) 
    c_edge_count_pos = defaultdict(int) 
    c_edge_count_chaos = defaultdict(int) 
    c_edge_count_total = defaultdict(int) 
    c_gene_presence_count = defaultdict(int) 
    c_gene_degree_count = defaultdict(int)
    c_background__degree = defaultdict(int)
    c_background__degree_prop = defaultdict(int)
    
    g2 = pickle.load(open(base_dir2, 'rb'))
    for node in g2.nodes():
        c_background__degree[node] = g2.degree(node)
    
    num_files = 0 # num patients 

    for f in glob.glob(net_dir + "/*.pickle"):
        num_files += 1 
        #print('[' + time.strftime("%H:%M:%S", time.localtime()) + "] Processing " + f + "....")
        g = pickle.load(open(f, 'rb'))

        for gene in g.nodes():
            gene_presence_count[gene] += 1
            gene_degree_count[gene] += g.degree(gene)
            
            c_gene_presence_count[gene] += 1
            c_gene_degree_count[gene] += g.degree(gene)
            if g.degree(gene) > g2.degree(gene):
                print(gene + ": g > g2")
            c_background__degree_prop[gene] += g.degree(gene)/g2.degree(gene)


        for gi, gj, data in g.edges(data=True):
            edge_count_total[gi] += 1
            edge_count_total[gj] += 1
            
            c_edge_count_total[gi] += 1
            c_edge_count_total[gj] += 1

            if data['chaos']:
                edge_count_chaos[gi] += 1
                edge_count_chaos[gj] += 1
                
                c_edge_count_chaos[gi] += 1
                c_edge_count_chaos[gj] += 1

            #edge_key = (gi, gj)
            if data['weight'] <= 0:
                #  negative edge
                edge_count_neg[gi] += 1
                edge_count_neg[gj] += 1
                
                c_edge_count_neg[gi] += 1
                c_edge_count_neg[gj] += 1
            else:
                #  positive edge
                edge_count_pos[gi] += 1
                edge_count_pos[gj] += 1
                
                c_edge_count_pos[gi] += 1
                c_edge_count_pos[gj] += 1
        
        metrics_per_cancer[can[0]] = {
            'edge_count_neg': c_edge_count_neg,
            'edge_count_pos': c_edge_count_pos,
            'edge_count_chaos': c_edge_count_chaos,
            'edge_count_total': c_edge_count_total,
            'gene_presence_count': c_gene_presence_count,
            'gene_degree_count': c_gene_degree_count,
            "original_background__degree": c_background__degree,
            "original_background__degree_norm": c_background__degree_prop, 
            'num_samples': num_files
            }


print("done!")

df_edge_count_neg = pd.DataFrame.from_dict(edge_count_neg, orient='index', columns=['edge_count_neg'])
df_edge_count_pos = pd.DataFrame.from_dict(edge_count_pos, orient='index', columns=['edge_count_pos'])
df_edge_count_chaos = pd.DataFrame.from_dict(edge_count_chaos, orient='index', columns=['edge_count_chaos'])
df_edge_count_total = pd.DataFrame.from_dict(edge_count_total, orient='index', columns=['edge_count_total'])
df_gene_presence_count = pd.DataFrame.from_dict(gene_presence_count, orient='index', columns=['gene_presence_count'])
df_gene_degree_count = pd.DataFrame.from_dict(gene_degree_count, orient='index', columns=['gene_degree_count'])
df_original_background__degree = pd.DataFrame.from_dict(gene_degree_count, orient='index', columns=['original_background__degree'])
df_presence_prop = pd.DataFrame()
df_degree_avg = pd.DataFrame()
df_degree_norm = pd.DataFrame()
df_degree_norm_prop = pd.DataFrame()


for cancer_type, metrics_dict in metrics_per_cancer.items():
    print(cancer_type)
    df_edge_count_neg[cancer_type + '_edge_count_neg'] = pd.Series(metrics_dict['edge_count_neg'])
    df_edge_count_pos[cancer_type + '_edge_count_pos'] = pd.Series(metrics_dict['edge_count_pos'])
    df_edge_count_chaos[cancer_type + '_edge_count_chaos'] = pd.Series(metrics_dict['edge_count_chaos'])
    df_edge_count_total[cancer_type + '_edge_count_total'] = pd.Series(metrics_dict['edge_count_total'])
    df_gene_presence_count[cancer_type + '_gene_presence_count'] = pd.Series(metrics_dict['gene_presence_count'])
    df_gene_degree_count[cancer_type + '_gene_degree_count'] = pd.Series(metrics_dict['gene_degree_count'])
    df_presence_prop[cancer_type + '_presence_prop'] = pd.Series({key: value / metrics_dict['num_samples'] if isinstance(value, (int, float)) else value 
               for key, value in metrics_dict['gene_presence_count'].items()}) # basically metrics_dict['gene_presence_count'] / metrics_dict['num_samples'] 
    df_degree_avg[cancer_type + '_degree_avg'] = pd.Series({key: value / metrics_dict['num_samples'] if isinstance(value, (int, float)) else value 
               for key, value in metrics_dict['gene_degree_count'].items()}) # basically metrics_dict['gene_degree_count'] / metrics_dict['num_samples']
    df_degree_norm[cancer_type + '_norm_degreed'] = pd.Series(metrics_dict['original_background__degree_norm'])
    df_degree_norm_prop[cancer_type + '_norm_degree_prop'] = pd.Series({key: metrics_dict['original_background__degree_norm'][key] / (metrics_dict['num_samples']) for key in metrics_dict['original_background__degree_norm']}) # basically metrics_dict['gene_degree_count'] / metrics_dict['original_background__degree']  

    

result_df = pd.concat([df_edge_count_neg, df_edge_count_pos, df_edge_count_chaos, df_edge_count_total,
                       df_gene_presence_count, df_gene_degree_count, df_presence_prop, df_degree_avg, df_degree_norm, df_degree_norm_prop], axis=1)
result_df.fillna(0, inplace=True)


rows = []
keys_list = metrics_per_cancer.keys()

for cancer_type in keys_list:
    print(cancer_type, end=",")
    column_name = cancer_type + '_gene_presence_count'
    column_degree = cancer_type + '_gene_degree_count'
    column_norm = cancer_type +'_norm_degree_prop'
    
    sorted_df = result_df.sort_values(by=[column_norm, column_name, column_degree], ascending=[False, False, False])
    
    top3_genes = sorted_df.head(10).index
    
    top3_gene_presence_counts = sorted_df.loc[top3_genes, column_name]
    top3_degreed = sorted_df.loc[top3_genes, column_degree]
    top3_norm = sorted_df.loc[top3_genes, column_norm]
    
    for rank, (gene, gene_presence_count, degree, norm) in enumerate(zip(top3_genes, top3_gene_presence_counts, top3_degreed, top3_norm), start=1):
        rows.append([gene, cancer_type, gene_presence_count, degree, norm, rank])

topgenes_df = pd.DataFrame(rows, columns=['gene', 'cancer_type', 'gene_presence_count', 'degree', 'norm', 'rank'])

assigned_genes = set()
row_list = []

for cancer_type in topgenes_df['cancer_type'].unique():
    outcancer_df = topgenes_df[topgenes_df['cancer_type'] == cancer_type]
    
    for _, row in outcancer_df.sort_values(by='rank').iterrows():
        gene = row['gene']
        
        if gene not in assigned_genes:
            assigned_genes.add(gene)
            row_list.append([gene, cancer_type])
            break  
            
top_genes = result_df.loc[result_df.index.isin(assigned_genes)]
top_genes.reset_index(inplace=True)
top_genes["index"] = top_genes["index"].astype(float)
top_genes_ent = pd.merge(top_genes, entrez_map, left_on='index', right_on='entrez', how='left')
top_genes_ent.set_index('symbol', inplace=True)
top_genes_ent = top_genes_ent.drop(columns=['uniprot','ensembl'])
top_genes_ent = top_genes_ent.drop_duplicates()

top_genes_ent.to_csv("/data/blue_heatmap.csv")


node_data = []
edge_data = []

for cancers in folders_tuples:
    print(cancers[0])
    ### fix
    path = f"/splitpea/{cancers[0]}/{cancers[1]}/output"

    file_list = [f for f in os.listdir(path) if f.endswith('.pickle')]
    for file_name in file_list:
        with open(os.path.join(path, file_name), 'rb') as f:
            net = pickle.load(f)
            node_degrees = net.degree()
            
            node_data.extend([
                {'node': node, 'degree': degree, 'sample': file_name[:28], 'cancer': cancers[0].split('_')[0]}
                for node, degree in node_degrees
            ])
            
            edge_data.extend([
                {'node1': edge[0], 'node2': edge[1], 'sample': file_name, 'cancer': cancers[0].split('_')[0]}
                for edge in net.edges()
            ])


symbol_to_entrez = dict(zip(entrez_map['symbol'], entrez_map['entrez']))
entrez_to_symbol_global = {value: key for key, value in symbol_to_entrez.items() if not pd.isna(value)}

results_cancer = []

file_path = '/ref/Cosmic_MutantCensus_v99_GRCh37.tsv.gz'
cosmic_genes_df = pd.read_csv(file_path, compression='gzip', sep='\t')

for cancers in folders_tuples:
    print('[' + time.strftime("%H:%M:%S", time.localtime()) + "] Loading networks... " + cancers[0])
    ### fix
    path = "/splitpea/" + cancers[0] + "/" + cancers[1] + "/output"


    file_list = os.listdir(path)
    for file_name in file_list:
        if file_name.endswith('.dat'):
            cosmic_genes_df_filter = cosmic_genes_df[cosmic_genes_df['SAMPLE_NAME'] == file_name[:15]]
            entrez_to_symbol = dict(zip(cosmic_genes_df_filter['entrez'], cosmic_genes_df_filter['GENE_SYMBOL']))

            entrez_set = set(cosmic_genes_df_filter['entrez'].dropna())
            #print(file_name)
            #net = pickle.load(open(path+ '/' + str(file_name), 'rb'))
            net_edges = pd.read_csv(path+ '/' + str(file_name), delim_whitespace=True)

            filtered_edges = net_edges[(net_edges['node1'].isin(entrez_set)) | (net_edges['node2'].isin(entrez_set))]

            for _, row in filtered_edges.iterrows():
                matched_entrez = row['node1'] if row['node1'] in entrez_set else row['node2']
                gene_symbol = entrez_to_symbol[matched_entrez]
                results_cancer.append({
                    'entrez': matched_entrez,
                    'GENE_SYMBOL': gene_symbol,
                    'node1': row['node1'],
                    'node2': row['node2'],
                    'weight': row['weight'],
                    'chaos': row['chaos'],
                    'sample': file_name[:28],
                    'cancer': cancers[0].split('_')[0]
                })

results_cancer_df = pd.DataFrame(results_cancer)
results_cancer_df['node1_symbol'] = results_cancer_df['node1'].map(entrez_to_symbol_global)
results_cancer_df['node2_symbol'] = results_cancer_df['node2'].map(entrez_to_symbol_global)
results_cancer_df['edge'] = list(zip(results_cancer_df['node1_symbol'], results_cancer_df['node2_symbol']))

mutated_driver_edges = results_cancer_df[["edge","sample","cancer"]]
mutated_driver_genes = results_cancer_df[["GENE_SYMBOL","sample","cancer"]]
driver_genes = set(mutated_driver_genes["GENE_SYMBOL"])


symbol_to_entrez = dict(zip(entrez_map['symbol'], entrez_map['entrez']))
entrez_to_symbol = {value: key for key, value in symbol_to_entrez.items() if not pd.isna(value)}

node_data_df = pd.DataFrame(node_data)
edge_data_df = pd.DataFrame(edge_data)

node_data_df['node'] = node_data_df['node'].astype(float)
node_data_df['GENE_SYMBOL'] = node_data_df['node'].map(entrez_to_symbol)
edge_data_df['node1'] = edge_data_df['node1'].astype(float)
edge_data_df['node2'] = edge_data_df['node2'].astype(float)
edge_data_df['node1_symbol'] = edge_data_df['node1'].map(entrez_to_symbol)
edge_data_df['node2_symbol'] = edge_data_df['node2'].map(entrez_to_symbol)
edge_data_df['edge'] = list(zip(edge_data_df['node1_symbol'], edge_data_df['node2_symbol']))
edge_data_df['sample'] = edge_data_df['sample'].str[:28]

print(node_data_df.head())
print(edge_data_df.head())

df_w_driver_info = edge_data_df.copy()
df_w_driver_info['node1_label'] = df_w_driver_info['node1_symbol'].apply(lambda x: 'driver' if x in driver_genes else 'other')
df_w_driver_info['node2_label'] = df_w_driver_info['node2_symbol'].apply(lambda x: 'driver' if x in driver_genes else 'other')

def process_sample(group):
    sample = group['sample'].iloc[0]
    cancer = group['cancer'].iloc[0]

    G = nx.Graph()
    edges = list(zip(group['node1'], group['node2']))
    G.add_edges_from(edges)

    node_labels = {}
    for _, row in group.iterrows():
        node_labels[row['node1']] = row['node1_label']
        node_labels[row['node2']] = row['node2_label']

    driver_labels = {'non-mutated driver', 'mutated driver', 'driver'}
    drivers = {node for node, label in node_labels.items() if label in driver_labels}

    if len(G) == 0:
        result = {'sample': sample, 'cancer': cancer,
                  '2_drivers': 0,
                  '1_driver': 0,
                  'no_drivers_connected_to_driver': 0,
                  'no_drivers_no_connection_to_driver': 0}
        return pd.Series(result)

    largest_cc = max(nx.connected_components(G), key=len)
    G_largest = G.subgraph(largest_cc).copy()

    node_has_driver_neighbor = {}
    for node in G_largest.nodes():
        node_has_driver_neighbor[node] = any(
            neighbor in drivers for neighbor in G_largest.neighbors(node)
        )

    counts = {
        '2_drivers': 0,
        '1_driver': 0,
        'no_drivers_connected_to_driver': 0,
        'no_drivers_no_connection_to_driver': 0
    }

    for node1, node2 in G_largest.edges():
        node1_is_driver = node_labels[node1] in driver_labels
        node2_is_driver = node_labels[node2] in driver_labels

        if node1_is_driver and node2_is_driver:
            counts['2_drivers'] += 1
        elif node1_is_driver or node2_is_driver:
            counts['1_driver'] += 1
        else:
            if node_has_driver_neighbor[node1] or node_has_driver_neighbor[node2]:
                counts['no_drivers_connected_to_driver'] += 1
            else:
                counts['no_drivers_no_connection_to_driver'] += 1

    result = {'sample': sample, 'cancer': cancer}
    result.update(counts)
    return pd.Series(result)


df_w_driver_info_diretvindirect = df_w_driver_info.groupby('sample').apply(process_sample).reset_index(drop=True)
columns_to_plot = ['2_drivers', '1_driver', 'no_drivers_connected_to_driver', 'no_drivers_no_connection_to_driver']
cancer_types_order = ["GBM", "LGG", "LAML", "DLBC", "TGCT", "ACC", "BLCA", "CESC", "UCEC", "UCS", "PAAD", "CHOL", "LIHC", "KICH", "KIRP", "KIRC",
                      "OV", "BRCA", "PRAD", "STAD", "THCA", "ESCA", "COAD", "READ", "SKCM", "MESO", "LUAD", "LUSC"]
res_melted = pd.melt(df_w_driver_info_diretvindirect, id_vars=['sample', 'cancer'], value_vars=columns_to_plot, 
                     var_name='driver_status', value_name='value')

res_melted['cancer'] = pd.Categorical(res_melted['cancer'], categories=cancer_types_order, ordered=True)
res_melted.to_csv("/data/driver_splitpea_info.csv")

mc3_mutated_genes = pd.read_csv('/home/jz94/splint/a_cleaned_code/external_data/mc3_mutation_data.csv')
driver_genes = set(mutated_driver_genes['GENE_SYMBOL'])

mutated_driver_genes["Participant"] = mutated_driver_genes["sample"].str.slice(0, 12)
node_data_df_copy = node_data_df.copy()
node_data_df_copy["Participant"] = node_data_df_copy["sample"].str.slice(0, 12)

merged_node_data_df = node_data_df_copy.merge(mutated_driver_genes,
                                  on=['Participant', 'GENE_SYMBOL'], how='left', indicator='driver_status')

mc3_mutated_genes['GENE_SYMBOL'] = mc3_mutated_genes['Hugo_Symbol']
merged_node_data_df = merged_node_data_df.merge(mc3_mutated_genes,
                                  on=['Participant', 'GENE_SYMBOL'], how='left', indicator='mc3_status')

def label_samples(row):
    if row['driver_status'] == 'both':
        return 'mutated driver'
    elif row['mc3_status'] == 'both' and row['GENE_SYMBOL'] not in driver_genes:
        return 'mutated non-driver'
    elif row['mc3_status'] == 'left_only' and row['GENE_SYMBOL'] in driver_genes:
        return 'non-mutated driver'
    else:
        return 'other'

merged_node_data_df['label'] = merged_node_data_df.apply(label_samples, axis=1)

driver_non_mutation_df = merged_node_data_df[['GENE_SYMBOL', 'Participant', 'degree', 'label']]
driver_non_mutation_df.to_csv('/data/driver_non_mutation_df.csv')