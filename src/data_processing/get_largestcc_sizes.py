#!/usr/bin/python

# takes a network + does some simple node property calculations 
# for the largest connected component

import networkx as nx

import pickle
import time

from pathlib import Path
import glob
import os

base_dir = str(Path(__file__).resolve().parent.parent)

def get_folders(directory):
    folders_list = []

    folders = [f for f in os.listdir(directory) if os.path.isdir(os.path.join(directory, f))]

    for folder in folders:
        subfolders = [f for f in os.listdir(os.path.join(directory, folder)) if os.path.isdir(os.path.join(directory, folder, f))]
        for subfolder in subfolders:
            folders_list.append((folder, subfolder))

    return folders_list


directory_path = "" #path to splitpea output directories
folders_tuples = get_folders(directory_path)
print(folders_tuples)

#net_dir = base_dir + "/IRIS/"
#cancers = ['BRCA', 'PAAD']

base_dir = ""

out = open("/results/tcga_splitpea_network.largest_ccs.txt", 'w')
out.write('\t'.join(["cancer", "sample", "direction", "orig_nodes", "orig_edges", "lcc_nodes", "lcc_edges"]) + '\n')

g = pickle.load(open("/splitpea/reference/human_ppi_ddi_bg.pickle", 'rb'))
g_largest_cc = g.subgraph(max(nx.connected_components(g), key=len)).copy()
out.write('\t'.join(["background",
                    "ppi_ddi_bg",
                    "all",
                    str(len(g)),
                    str(g.number_of_edges()),
                    str(len(g_largest_cc)),
                    str(g_largest_cc.number_of_edges())]) + '\n')

for c in folders_tuples:
    print('[' + time.strftime("%H:%M:%S", time.localtime()) + "] Loading " + c[0] + " networks...")

    net_dir_ = base_dir + c[0] + '/' + c[1] + '/output/'
    for f in os.listdir(net_dir_):
        if 'edges.pickle' in f:
            file_pref = f.split('/')[-1].split('.')[0]
            g = pickle.load(open(net_dir_+ f, 'rb'))
            
            g_largest_cc = g.subgraph(max(nx.connected_components(g), key=len)).copy()

            out.write('\t'.join([c[0],
                                file_pref,
                                "all",
                                str(len(g)),
                                str(g.number_of_edges()),
                                str(len(g_largest_cc)),
                                str(g_largest_cc.number_of_edges())]) + '\n')
        
            g_largest_cc_neg = g_largest_cc.copy()
            g_largest_cc_pos = g_largest_cc.copy()
            g_largest_cc_chaos = g_largest_cc.copy()
            
            for gi, gj in g_largest_cc.edges():
                if g_largest_cc[gi][gj]['chaos']:
                    g_largest_cc_neg.remove_edge(gi, gj)
                    g_largest_cc_pos.remove_edge(gi, gj)
                    continue
                else:
                    g_largest_cc_chaos.remove_edge(gi, gj)

                if g_largest_cc_pos[gi][gj]['weight'] <= 0:
                    g_largest_cc_pos.remove_edge(gi, gj)
                else:
                    g_largest_cc_pos[gi][gj]['distance'] = 1 - g_largest_cc_pos[gi][gj]['weight']
                    g_largest_cc_neg.remove_edge(gi, gj)
            
            g_largest_cc_neg.remove_nodes_from([n for (n, deg) in g_largest_cc_neg.degree() if deg == 0]) 
            g_largest_cc_pos.remove_nodes_from([n for (n, deg) in g_largest_cc_pos.degree() if deg == 0])
            g_largest_cc_chaos.remove_nodes_from([n for (n, deg) in g_largest_cc_chaos.degree() if deg == 0]) 

            out.write('\t'.join([c[0],
                                file_pref,
                                "negative",
                                str(len(g_largest_cc)),
                                str(g_largest_cc.number_of_edges()),
                                str(len(g_largest_cc_neg)),
                                str(g_largest_cc_neg.number_of_edges())]) + '\n')
        
            out.write('\t'.join([c[0],
                                file_pref,
                                "positive",
                                str(len(g_largest_cc)),
                                str(g_largest_cc.number_of_edges()),
                                str(len(g_largest_cc_pos)),
                                str(g_largest_cc_pos.number_of_edges())]) + '\n')

            out.write('\t'.join([c[0],
                                file_pref,
                                "chaos",
                                str(len(g_largest_cc)),
                                str(g_largest_cc.number_of_edges()),
                                str(len(g_largest_cc_chaos)),
                                str(g_largest_cc_chaos.number_of_edges())]) + '\n')
        
out.close()