#!/usr/bin/python
import pandas as pd
import os  
import genopyc as gp
import networkx as nx


def create_inputs():

    # Load the data and parse the interactome
    interactome = pd.read_csv('../data/HIPPIE-current.mitab.txt',sep='\t',low_memory=False)
    interactome = interactome[interactome.columns[:2]]
    interactome[interactome.columns[0]] = interactome[interactome.columns[0]].apply(lambda x: x.replace('entrez gene:',''))
    interactome[interactome.columns[1]] = interactome[interactome.columns[1]].apply(lambda x: x.replace('entrez gene:',''))
    interactome = interactome[~interactome['ID Interactor A'].str.contains('-')]
    interactome = interactome[~interactome['ID Interactor B'].str.contains('-')]
    interactome = interactome.astype(int)
    interactome['ID_A_symbol'] = gp.geneId_mapping(interactome['ID Interactor A'].tolist(),'entrez','symbol')
    interactome['ID_B_symbol'] = gp.geneId_mapping(interactome['ID Interactor B'].tolist(),'entrez','symbol')
    interactome.dropna(inplace=True)
    interactome_net = nx.from_pandas_edgelist(interactome[['ID_A_symbol','ID_B_symbol']],source='ID_A_symbol',target='ID_B_symbol')
    largest_cc = max(nx.connected_components(interactome_net), key=len)
    interactome_cc = nx.subgraph(interactome_net,largest_cc)
    interactome_pandas = nx.to_pandas_edgelist(interactome_cc)
    interactome_pandas.insert(1,'score',1)
    il4_edges = interactome[(interactome.ID_A_symbol == 'IL4') | (interactome.ID_B_symbol == 'IL4')].iloc[:,2:].values
    il10_edges = interactome[(interactome.ID_A_symbol == 'IL10') | (interactome.ID_B_symbol == 'IL10')].iloc[:,2:].values


    Dipper = pd.read_csv('../data/dipper_protein_expression.csv')
    Dipper=Dipper.iloc[:,2:]
    Dipper.columns=[x.replace(".", "_") for x in Dipper.columns]
    columnyoungfilter=[column for column in Dipper.columns.tolist() if ('young_NP' in column) ]
    columnoldfilter=[column for column in Dipper.columns.tolist() if ('old_NP' in column) ]

    columnyoungfilter=columnyoungfilter+list(Dipper.columns[:3])
    columnoldfilter = columnoldfilter+list(Dipper.columns[:3])

    dipper_young = Dipper[columnyoungfilter]
    dipper_old = Dipper[columnoldfilter]
    dipper_young = dipper_young.dropna(subset=['L3_4_young_NP','L4_5_young_NP','L5_S1_young_NP'],how='all')
    dipper_old = dipper_old.dropna(subset=['L3_4_old_NP','L4_5_old_NP','L5_S1_old_NP'],how='all')

    gene_symbol_list_seeds=['IL6', 'IL6ST','IL16', 'IL1B', 'IL17D', 'IL18', 'IL20', 'TNF', 'LIF', 'OSM', 'CCL2', 'CCL3', 'CCL5', 'CCL7', 'CXCL1', 'CXCL2', 'CXCL3', 'CXCL8', 'CX3CL1', 'IGF1', 'GDF5', 'TGFB1', 'GDF6', 'PMP2', 'VEGFA']
    #CCL4 was removed due to the absence of data in the interactome
    
    # get the genes of dipper young and old
    genes_dipper_young = dipper_young['Gene_names'].tolist()
    genes_dipper_old = dipper_old['Gene_names'].tolist()

    # The seeds have been found in the lab to be expressed in the disc so we will force their presence in the final interactomes
    missing_seeds_young = [gene for gene in gene_symbol_list_seeds if gene not in genes_dipper_young]
    missing_seeds_old = [gene for gene in gene_symbol_list_seeds if gene not in genes_dipper_old]
    genes_interactome_expressed_young = dipper_young.Gene_names.tolist() + missing_seeds_young
    genes_interactome_expressed_old = dipper_old.Gene_names.tolist() + missing_seeds_old
    genes_interactome_expressed_young = sum([a.split(';') for a in genes_interactome_expressed_young],[])
    genes_interactome_expressed_old = sum([a.split(';') for a in genes_interactome_expressed_old],[])

    # We will force also in the interactome to be present the genes that interact with Il4 and IL10 that are expressed 
    filtered_connections_il10 = [edge for edge in il10_edges if ('IGHG1' in edge) or ('A2M' in edge)]
    filtered_connections_il4 = [edge for edge in il4_edges if 'A2M' in edge]
    filtered_connections_il10_df = pd.DataFrame([[a,1,b] for [a,b] in filtered_connections_il10],columns=['source','score','target'])
    filtered_connections_il4_df = pd.DataFrame([[a,1,b] for [a,b] in filtered_connections_il4],columns=['source','score','target'])

    # Create the filtered interactomes
    interactome_young = interactome_pandas[(interactome_pandas.source.isin(genes_interactome_expressed_young)) | (interactome_pandas.target.isin(genes_interactome_expressed_young))]
    interactome_old = interactome_pandas[(interactome_pandas.source.isin(genes_interactome_expressed_old)) | (interactome_pandas.target.isin(genes_interactome_expressed_old))]
    interactome_young_IL4 = pd.concat([interactome_young,filtered_connections_il4_df])
    interactome_old_IL4 = pd.concat([interactome_old,filtered_connections_il4_df])
    interactome_young_IL10 = pd.concat([interactome_young,filtered_connections_il10_df])
    interactome_old_IL10 = pd.concat([interactome_old,filtered_connections_il10_df])
    

    # Make sure the filtered young and old interactomes are connected
    interactome_young_net = nx.from_pandas_edgelist(interactome_young)
    interactome_old_net = nx.from_pandas_edgelist(interactome_old)
    interactome_young_net_IL4 = nx.from_pandas_edgelist(interactome_young_IL4)
    interactome_old_net_IL4 = nx.from_pandas_edgelist(interactome_old_IL4)
    interactome_young_net_IL10 = nx.from_pandas_edgelist(interactome_young_IL10)
    interactome_old_net_IL10 = nx.from_pandas_edgelist(interactome_old_IL10)

    # Create node file
    nodes_young = list(set(interactome_young.source.tolist() + interactome_young.target.tolist()))
    nodes_old = list(set(interactome_old.source.tolist() + interactome_old.target.tolist()))
    nodes_old_IL4 = list(set(interactome_old_IL4.source.tolist() + interactome_old_IL4.target.tolist()))
    nodes_old_IL10 = list(set(interactome_old_IL10.source.tolist() + interactome_old_IL10.target.tolist()))
    nodes_young_IL4 = list(set(interactome_young_IL4.source.tolist() + interactome_young_IL4.target.tolist()))
    nodes_young_IL10 = list(set(interactome_young_IL10.source.tolist() + interactome_young_IL10.target.tolist()))

    # Save the interactomes
    def get_var_name(variable):
        for name, value in globals().items():
            if value is variable:
                return name

    #check if the interactomes are all connected
    interactomes = [interactome_old_net_IL10 ,interactome_young_net_IL10,interactome_old_net_IL4,interactome_young_net_IL4,interactome_old_net,interactome_young_net]

    for inte in interactomes:
        int_name = get_var_name(inte)
        largest_cc = max(nx.connected_components(inte), key=len)
        if len(largest_cc)==len(inte.nodes()):
            pass
        else:
            print(f"The there is discrepancies between lcc and interactome for interactome {int_name}")


    # Save the interactomes
    interactomes = [interactome_old_IL10 ,interactome_young_IL10, interactome_old_IL4,interactome_young_IL4,interactome_old,interactome_young]
    for inte in interactomes:
        inte_name =  get_var_name(inte)
        inte.to_csv(f'../inputs/Dipper/{inte_name}.txt',index=False,header=False,sep = '\t')

    # Create seeds files
    seeds_young = list(zip(nodes_young, [1 if gene in gene_symbol_list_seeds else 0 for gene in nodes_young]))
    # write the inputs
    with open(f'../inputs/Dipper/seeds_young.txt','w') as f:
        for s in seeds_young:
            f.write(str(s[0]) + '\t' + str(s[1]) + '\n')

    seeds_young_IL4 = list(zip(nodes_young_IL4, [1 if gene in gene_symbol_list_seeds +['IL4'] else 0 for gene in nodes_young_IL4]))
    # write the inputs
    with open(f'../inputs/Dipper/seeds_young_IL4.txt','w') as f:
        for s in seeds_young_IL4:
            f.write(str(s[0]) + '\t' + str(s[1]) + '\n')
            
    seeds_young_IL10 = list(zip(nodes_young_IL10, [1 if gene in gene_symbol_list_seeds +['IL10'] else 0 for gene in nodes_young_IL10]))
    # write the inputs
    with open(f'../inputs/Dipper/seeds_young_IL10.txt','w') as f:
        for s in seeds_young_IL10:
            f.write(str(s[0]) + '\t' + str(s[1]) + '\n')

    seeds_old= list(zip(nodes_old, [1 if gene in gene_symbol_list_seeds else 0 for gene in nodes_old]))
    # write the inputs
    with open(f'../inputs/Dipper/seeds_old.txt','w') as f:
        for s in seeds_old:
            f.write(str(s[0]) + '\t' + str(s[1]) + '\n')

    seeds_old_IL4 = list(zip(nodes_old_IL4, [1 if gene in gene_symbol_list_seeds +['IL4'] else 0 for gene in nodes_old_IL4]))
    # write the inputs
    with open(f'../inputs/Dipper/seeds_old_IL4.txt','w') as f:
        for s in seeds_old_IL4:
            f.write(str(s[0]) + '\t' + str(s[1]) + '\n')
            
    seeds_old_IL10 = list(zip(nodes_old_IL10, [1 if gene in gene_symbol_list_seeds +['IL10'] else 0 for gene in nodes_old_IL10]))
    # write the inputs
    with open(f'../inputs/Dipper/seeds_old_IL10.txt','w') as f:
        for s in seeds_old_IL10:
            f.write(str(s[0]) + '\t' + str(s[1]) + '\n')
            
    # Save the interactomes
    def get_var_name(variable):
        for name, value in locals().items():
            if value is variable:
                return name

    #check if the interactomes are all connected
    interactomes = [interactome_old_net_IL10 ,interactome_young_net_IL10,interactome_old_net_IL4,interactome_young_net_IL4,interactome_old_net,interactome_young_net]

    for inte in interactomes:
        int_name = get_var_name(inte)
        largest_cc = max(nx.connected_components(inte), key=len)
        if len(largest_cc)==len(inte.nodes()):
            pass
        else:
            print(f"The there is discrepancies between lcc and interactome for interactome {int_name}")



    interactomes = [interactome_old_IL10 ,interactome_young_IL10, interactome_old_IL4,interactome_young_IL4,interactome_old,interactome_young]
    for inte in interactomes:
        inte_name =  get_var_name(inte)
        inte.to_csv(f'../inputs/Dipper/{inte_name}_txt',index=False,header=False,sep = '\t')
        


def run_guild():
    experiments = ['All_young','All_old','IL4_young','IL4_old','IL10_young','IL10_old']
    for i,exp in enumerate(experiments):
        # The interactome have to be composed from only the genes experessed either in old and young dipper
        if i == 0:
            print(f"Performing experiment {exp}")
            # running guild
            os.system(f'../guild/guild_x64 -s s -n ../inputs/Dipper/seeds_young.txt -e ../inputs/Dipper/interactome_young.txt -o ../outputs/Dipper/{exp}_Dipper_out.txt -i 2 -r 3')
        elif i == 1:
            print(f"Performing experiment {exp}")
            # run guild
            os.system(f'../guild/guild_x64 -s s -n ../inputs/Dipper/seeds_old.txt -e ../inputs/Dipper/interactome_old.txt -o ../outputs/Dipper/{exp}_Dipper_out.txt -i 2 -r 3')
            
        elif i == 2:
            print(f"Performing experiment {exp}")
            # run guild
            os.system(f'../guild/guild_x64 -s s -n ../inputs/Dipper/seeds_young_IL4.txt -e ../inputs/Dipper/interactome_young_IL4.txt -o ../outputs/Dipper/{exp}_Dipper_out.txt -i 2 -r 3')
        elif i == 3:
            print(f"Performing experiment {exp}")
            # run guild
            os.system(f'../guild/guild_x64 -s s -n ../inputs/Dipper/seeds_old_IL4.txt -e ../inputs/Dipper/interactome_old_IL4.txt -o ../outputs/Dipper/{exp}_Dipper_out.txt -i 2 -r 3')
        elif i == 4:
            print(f"Performing experiment {exp}")
            # run guild
            os.system(f'../guild/guild_x64 -s s -n ../inputs/Dipper/seeds_young_IL10.txt -e ../inputs/Dipper/interactome_young_IL10.txt -o ../outputs/Dipper/{exp}_Dipper_out.txt -i 2 -r 3')
        elif i == 5:
            print(f"Performing experiment {exp}")
            # run guild
            os.system(f'../guild/guild_x64 -s s -n ../inputs/Dipper/seeds_old_IL10.txt -e ../inputs/Dipper/interactome_old_IL10.txt -o ../outputs/Dipper/{exp}_Dipper_out.txt -i 2 -r 3')


def Main():
    create_inputs()
    run_guild()
    
if __name__ == '__main__':
    Main()
    
                
           
       
                

    





    
    
    
    
    
