#!/usr/bin/python
import pandas as pd
import os  
import genopyc as gp
import networkx as nx


def minmax_scale_list(data, min_value=0, max_value=1):
    """
    Min-max scale a list of numbers to a specified range.

    Parameters:
        data (list): The list of numbers to be scaled.
        min_value (float): The minimum value of the desired range (default: 0).
        max_value (float): The maximum value of the desired range (default: 1).

    Returns:
        list: The scaled list of numbers.
    """
    min_data = min(data)
    max_data = max(data)
    scaled_data = [(x - min_data) / (max_data - min_data) * (max_value - min_value) + min_value for x in data]
    return scaled_data

def load_data():
    secretome_deg = pd.read_excel('../data/CORRECT_Secretome profile human degenerated NPC.xlsx')
    secretome_healthy = pd.read_excel('../data/Secretome profile human trauma NPC.xlsx')
    cols_secretome = secretome_deg.iloc[3]
    secretome_healthy.columns = cols_secretome
    secretome_deg.columns = cols_secretome
    secretome_healthy = secretome_healthy[4:]
    secretome_deg = secretome_deg[4:]
    secretome_deg = secretome_deg.iloc[58:]
    secretome_healthy = secretome_healthy.iloc[:9]
    secretome_deg['Analyte:'] = secretome_deg['Analyte:'].apply(lambda x: x.replace('+',''))
    secretome_healthy['Analyte:'] = secretome_healthy['Analyte:'].apply(lambda x: x.replace('+',''))
    secretome_deg.set_index('Analyte:',inplace=True)    
    secretome_healthy.set_index('Analyte:',inplace=True)
    secretome_deg = secretome_deg.drop(['CCL4', 'IL8'],axis=1)
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
    return secretome_deg,secretome_healthy,interactome

def parse_network(interactome):
    interactome_net = nx.from_pandas_edgelist(interactome[['ID_A_symbol','ID_B_symbol']],source='ID_A_symbol',target='ID_B_symbol')
    largest_cc = max(nx.connected_components(interactome_net), key=len)
    interactome_cc = nx.subgraph(interactome_net,largest_cc)
    interactome_pandas = nx.to_pandas_edgelist(interactome_cc)
    interactome_pandas.insert(1,'score',1)
    return interactome_pandas

def write_guild_inputs_deg(secretome_deg):
    
    experiments = ['-IL10','IL1IL4','IL1IL10']

    for i,exp in enumerate(experiments):
        df = secretome_deg[secretome_deg.index.str.contains(exp)]
        
        means = [df[col].mean() for col in df.columns.tolist()]
        scaled_means = minmax_scale_list(means)
        
        # Make a dictionary of cytokynes and values
        dictionary_of_means = dict(zip(df.columns.tolist(),scaled_means))
        
        # Force manually depending on the experiment
        if i == 0:
            dictionary_of_means['IL4'] = 1
        elif i == 1:
            dictionary_of_means['IL4'] = 1
            dictionary_of_means['IL1'] = 1
        elif i == 2:
            dictionary_of_means['IL10'] = 1
            dictionary_of_means['IL1'] = 1

        # Write the files
        with open('../inputs/secretome/degenerated_'+exp+'_nodes.txt','w') as f:
            for gene,score in dictionary_of_means.items():
                f.write(gene + ' ' + str(score) + ' ' + '\n')
                
def write_guild_inputs_healthy(secretome_healthy):
    experiments = ['Ctrl','IL4','IL1']
    
    for i,exp in enumerate(experiments):
        
        df = secretome_healthy[secretome_healthy.index.str.contains(exp)]
        means = [df[col].mean() for col in df.columns.tolist()]
        scaled_means = minmax_scale_list(means)
        dictionary_of_means = dict(zip(df.columns.tolist(),scaled_means))
        if i == 1:
            dictionary_of_means['IL4'] = 1
        elif i == 2:
            dictionary_of_means['IL1'] = 1
            
            # Write the files
        with open('../inputs/secretome/healthy_'+exp+'_nodes.txt','w') as f:
            for gene,score in dictionary_of_means.items():
                f.write(gene + ' ' + str(score) + ' ' + '\n')      


def run_guild():
    interactome_path = '../inputs/'
    nodes_path = '../inputs/secretome/'
    output_path = '../outputs/secretome/'
    
    for file in os.listdir('../inputs/secretome/'):
        exp_name = file.replace('_nodes.txt','')
        os.system('chmod 777 ../inputs/secretome/')
        os.system(f'../guild/guild_x64 -s s -n {nodes_path+file} -e {interactome_path}interactome.sif -o {output_path+exp_name} -i 2 -r 3')          

def Main():
    print('=========================== LOADING DATA ===========================')
    secretome_deg, secretome_healthy, interactome = load_data()
    print('=========================== PARSING INTERACTOME ===========================')
    interactome_pandas = parse_network(interactome)
    print('=========================== WRITING INPUT FILES ===========================')
    interactome_pandas.to_csv('../inputs/interactome.sif',index=False,header=False)
    write_guild_inputs_deg(secretome_deg)
    write_guild_inputs_healthy(secretome_healthy)
    print('=========================== RUNNING GUILD ===========================')

    run_guild()
    
    
if __name__ == '__main__':
    Main()
    

        
        

    
    
    
