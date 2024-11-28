import networkx as nx
import pandas as pd
import os
import pickle
import graphlot as gp

def Main():
    dipper_interactome_path='../inputs/Dipper/'
    secretome_interactome_path='../inputs/subnet_enriched.sif'
    explants_interactome_path='../inputs/subnet_enriched.sif'
    seeds_secretome_path = '../inputs/secretome/'
    seeds_secretome_path = '../inputs/explants/'
    table = pd.DataFrame(columns=["N_of_nodes","Diameter","Center"])
    res = []
    for net_name in os.listdir(dipper_interactome_path):
        if 'interactome' in net_name:
            sif_file = pd.read_csv(dipper_interactome_path+net_name,sep='\t',names=['source','score','target'])
            net = nx.from_pandas_edgelist(sif_file)
            connected_components = list(nx.connected_components(net))
            largest_component = max(connected_components, key=len)
            net_conn = net.subgraph(largest_component)
            n_of_nodes = net_conn.number_of_nodes()
            diameter = nx.diameter(net_conn)
            center = nx.center(net_conn)
        res.append([net_name,n_of_nodes,diameter,center])
    with open("../outputs/interactome_stats.pickle",'rb') as f:
        pickle.dump(res,f)


if __name__ == '__main__':
    Main()