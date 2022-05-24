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
if not os.path.isdir('./Data/Clustes'): os.mkdir('./Data/Clusters')

#%% Loading the necessary data

with open('./Data/Metadata/metadata.json') as file:
    metadata = json.load(file)

#%% All clusterings

# Done on HPC: Cluster.sh

# with open('./Data/Clusteres/Cluster_INFO.txt', 'a') as info_file:
#     for identity in ['0.9995', '0.9985', '0.998','0.9975', '0.995', '0.99']:
#         cluster_count = 0
#         seq_counts = []
#         with open('./Data/Clustered/SARSCoV2_'+identity[2:]+'.fasta.clstr', 'r') as file:
#             line_list = file.readlines()
#             for count in range(0, len(line_list)):
#                 line = line_list[count]
#                 if line.startswith('>Cluster'):
#                     cluster_count += 1
#                     if cluster_count >= 1:
#                         seq_counts.append(line_list[count-1].split()[0])
#         seq_counts.append(line_list[-1].split()[0])
#         seq_counts = list(map(int, seq_counts))
#         info_file.write('Identity: {} \n\t\tclusters: {}\n\t\tmax_seq: {}\n\t\tmin_seq: {}\n'.format(identity, cluster_count, max(seq_counts), min(seq_counts)))
#         print('Identity: {} \n\t\tclusters: {}\n\t\tmax_seq: {}\n\t\tmin_seq: {}\n'.format(identity, cluster_count, max(seq_counts), min(seq_counts)))

# Identity: 0.9995 
# 		clusters: 14684
# 		max_seq: 58819
# 		min_seq: 0

# Identity: 0.9985 
# 		clusters: 301
# 		max_seq: 148697
# 		min_seq: 0

# Identity: 0.998 
# 		clusters: 185
# 		max_seq: 243443
# 		min_seq: 0

# Identity: 0.9975 
# 		clusters: 143
# 		max_seq: 246232
# 		min_seq: 0

# Identity: 0.995 
# 		clusters: 81
# 		max_seq: 248315
# 		min_seq: 0

# Identity: 0.99 
# 		clusters: 57
# 		max_seq: 249218
# 		min_seq: 0

# => 0.9985 is probably best

#%% Extracting virus name from sequence description

cluster_metadata = {}

for seq_record in SeqIO.parse('./Data/Clusters/SARSCoV2_9985.fasta', 'fasta'):
    ID = seq_record.id
    cluster_metadata[ID] = metadata[ID]

for ID in cluster_metadata:
    des = cluster_metadata[ID]['Description']
    if 'assembly' in des:
        name = ID + ' genome assembly'
    else: 
        des = des.split()
        for word in des:
            if any([i in word for i in ['SARS-CoV-2/', 'SARS-Cov-2/', 'hCoV-19/', 'BetaCoV/']]):
                   name = word.strip(',')
    cluster_metadata[ID]['Name'] = name

with open('./Data/Metadata/cluster_metadata.json', 'w') as f:
    json.dump(cluster_metadata, f)

#%% Making new fasta file with aligned sequence and virus name

with open('./Data/Tree/SARSCoV2.fasta', 'w') as f:
    for seq_record in SeqIO.parse('./Data/Alignment/SARSCoV2.fasta', 'fasta'):
        ID = seq_record.id
        seq = str(seq_record.seq)
        if ID in cluster_metadata:
            name = cluster_metadata[ID]['Name']
            f.write('>' + name + '\n' + seq + '\n')
    
#%% Creating a tree via FastTree

input_file = './Data/Tree/SARSCoV2.fasta'
output_file = './Data/Tree/SARSCoV2.tree'
cmd = 'fasttree -nt -quote '+input_file+' > '+output_file
subprocess.run(cmd, shell = True)
