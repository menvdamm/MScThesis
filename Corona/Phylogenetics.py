#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 15 13:32:43 2021

@author: mennovandamme
"""

### Clustering & Phylogentic tree construction

#%% Dependencies

# conda install -c bioconda cd-hit
# conda install -c bioconda fasttree
import os, subprocess, json
from Bio import SeqIO

#%% Directories

if not os.path.isdir('./Data/Tree'): os.mkdir('./Data/Tree')
if not os.path.isdir('./Data/Clusters'): os.mkdir('./Data/Clusters')

#%% Loading the necessary data

with open('./Data/Metadata/metadata.json') as file:
    metadata = json.load(file)

#%% All clusterings

# Done on HPC: Cluster.sh

with open('./Data/Clusters/Cluster_INFO.txt', 'w') as info_file:
    for identity in ['0.9995', '0.999', '0.9985', '0.998','0.9975', '0.995', '0.99']:
        cluster_count = 0
        with open('./Data/Clusters/SARSCoV2_'+identity[2:]+'.fasta.clstr', 'r') as f:
            for line in f:
                if line.startswith('>Cluster'):
                    cluster_count += 1
        info_file.write('Identity: {} \n\tclusters: {}\n'.format(identity, cluster_count))

# Identity: 0.9995 
# 		clusters: 14684
# Identity: 0.999 
# 	clusters: 1117
# Identity: 0.9985 
# 		clusters: 301
# Identity: 0.998 
# 		clusters: 185
# Identity: 0.9975 
# 		clusters: 143
# Identity: 0.995 
# 		clusters: 81
# Identity: 0.99 
# 		clusters: 57

#%% Extracting virus name from sequence description

cluster_metadata = {}

for seq_record in SeqIO.parse('./Data/Clusters/SARSCoV2_9985.fasta', 'fasta'):
    ID = seq_record.id
    cluster_metadata[ID] = metadata[ID]


for ID in cluster_metadata:
    des = cluster_metadata[ID]['Description']
    if any([i in des for i in ['SARS-CoV-2/', 'SARS-Cov-2/', 'hCoV-19/', 'BetaCoV/']]):
        des = des.split()
        for word in des:
            if any([i in word for i in ['SARS-CoV-2/', 'SARS-Cov-2/', 'hCoV-19/', 'BetaCoV/']]):
                name = word.strip(',')
                splitname = name.split('/')  
                if splitname[1] == 'human':
                    name = '/'.join(name.split('/')[2:])
                elif splitname[0] == 'hCoV-19':
                    name = '/'.join(name.split('/')[1:])
    else:
        name = ID
    name += '(' + cluster_metadata[ID]['Pangolin_class'] + ')'
    cluster_metadata[ID]['Name'] = name

with open('./Data/Metadata/cluster_metadata.json', 'w') as f:
    json.dump(cluster_metadata, f)

#%% Making new fasta file with aligned sequence and virus name

with open('./Data/Tree/SARSCoV2.fasta', 'w') as f:
    for seq_record in SeqIO.parse('./Data/Alignment/SARSCoV2.fasta', 'fasta'):
        ID = seq_record.id
        if ID in cluster_metadata:
            seq = str(seq_record.seq)
            name = cluster_metadata[ID]['Name']
            f.write('>' + name + '\n' + seq + '\n')
    
#%% Creating a tree via FastTree

input_file = './Data/Tree/SARSCoV2.fasta'
output_file = './Data/Tree/SARSCoV2.tree'
cmd = 'fasttree -nt -quote '+input_file+' > '+output_file
subprocess.run(cmd, shell = True)
