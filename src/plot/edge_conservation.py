import pandas as pd
import matplotlib.pyplot as plt 


df_analysis = pd.read_csv('/results/edge_stats.csv') # from consensus_analysis.py

thresholds = [0, 0.25, .5, .75, 1]

counts_per_edge_pos = pd.DataFrame()
counts_per_edge_neg = pd.DataFrame()

for thresh in thresholds:
    print(thresh,end = "")
    if thresh != 0:
        df_thresh = df_analysis[df_analysis['Positive_Rate'] >= thresh]
    else:
        df_thresh = df_analysis[df_analysis['Positive_Rate'] > thresh]
    
    counts = df_thresh.groupby('Edge')['cancer'].nunique().reset_index()
    counts.rename(columns={'cancer': 'Cancer_Count'}, inplace=True)
    counts['Threshold'] = thresh  
    print(".", end="")
    
    counts_per_edge_pos = pd.concat([counts_per_edge_pos, counts], ignore_index=True)
    
    if thresh != 0:
        df_thresh = df_analysis[df_analysis['Negative_Rate'] >= thresh]
    else:
        df_thresh = df_analysis[df_analysis['Negative_Rate'] > thresh]
    
    counts = df_thresh.groupby('Edge')['cancer'].nunique().reset_index()
    counts.rename(columns={'cancer': 'Cancer_Count'}, inplace=True)
    counts['Threshold'] = thresh  
    print(".")
    
    counts_per_edge_neg = pd.concat([counts_per_edge_neg, counts], ignore_index=True)


freq_distribution = counts_per_edge_pos.groupby(['Threshold', 'Cancer_Count']).size().reset_index(name='Edge_Count')

plt.figure(figsize=(8, 6))
for thresh in thresholds:
    subset = freq_distribution[freq_distribution['Threshold'] == thresh]
    plt.plot(subset['Cancer_Count'], subset['Edge_Count'], marker='o', label=f'{thresh}')

plt.xlabel('Number of Cancers Tissue Pairs per Edge')
plt.ylabel('Number of Edges')
plt.legend(title='Positive Threshold')
plt.grid(True)
plt.savefig('/home/jz94/splint/a_cleaned_code/plots/fig4/pos_cons_cons.pdf', bbox_inches='tight')
plt.show()

freq_distribution = counts_per_edge_neg.groupby(['Threshold', 'Cancer_Count']).size().reset_index(name='Edge_Count')
plt.figure(figsize=(8, 6))
for thresh in thresholds:
    subset = freq_distribution[freq_distribution['Threshold'] == thresh]
    plt.plot(subset['Cancer_Count'], subset['Edge_Count'], marker='o', label=f'{thresh}')
plt.xlabel('Number of Cancers Tissue Pairs per Edge')
plt.ylabel('Number of Edges')
plt.legend(title='Negative Threshold')
plt.grid(True)
plt.savefig('/home/jz94/splint/a_cleaned_code/plots/fig4/neg_cons_cons.pdf', bbox_inches='tight')
plt.show()