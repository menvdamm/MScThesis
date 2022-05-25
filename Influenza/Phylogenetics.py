#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 15 13:32:43 2021

@author: mennovandamme
"""

#%% Dependencies

# conda install -c bioconda cd-hit
# conda install -c bioconda fasttree
import os, subprocess, json, re
from Bio import SeqIO
import pandas as pd

#%% Directories

if not os.path.isdir('./Data/Clusters'): os.mkdir('./Data/Clusters')
if not os.path.isdir('./Data/Tree'): os.mkdir('./Data/Tree')

#%% Loading the necessary data

with open('./Data/Metadata/metadata.json') as file:
    metadata = json.load(file)
    
#%% All clusterings

cluster_info = pd.DataFrame({'Type': ['A']*8 + ['B']*8,
                             'Segment': list('12345678')*2}) 

for identity in ['0.995', '0.99', '0.98', '0.97', '0.96', '0.95']:
    cluster_counts = []
    for Type in ['A', 'B']:
        for Segment in list('12345678'):     
            input_file = './Data/Clean/'+Type+'_'+Segment+'.fasta'
            output_file = './Data/Clusters/'+Type+'_'+Segment+'.fasta'
            cmd = 'cd-hit-est -i '+input_file+' -o '+output_file+' -c '+identity+' -n 10 -T 0'
            subprocess.run(cmd, shell = True)
            with open('./Data/Clusters/'+Type+'_'+Segment+'.fasta.clstr', 'r') as f:
                cluster_count = 0
                for line in f:
                    if line.startswith('>Cluster'):
                        cluster_count += 1
            cluster_counts.append(cluster_count)
    cluster_info[str(float(identity)*100)] = cluster_counts
                
cluster_info.to_csv('./Data/Clusters/cluster_info.csv', index = False)                   

#%% Clustering

for Type in ['A']:
    for Segment in list('12345678'):
        identity = '0.99'
        input_file = './Data/Clean/'+Type+'_'+Segment+'.fasta'
        output_file = './Data/Clusters/'+Type+'_'+Segment+'.fasta'
        cmd = 'cd-hit-est -i '+input_file+' -o '+output_file+' -c '+identity+' -n 10 -T 0'
        subprocess.run(cmd, shell = True)
        
for Type in ['B']:
    for Segment in list('12345678'):
        identity = '0.995'
        input_file = './Data/Clean/'+Type+'_'+Segment+'.fasta'
        output_file = './Data/Clusters/'+Type+'_'+Segment+'.fasta'
        cmd = 'cd-hit-est -i '+input_file+' -o '+output_file+' -c '+identity+' -n 10 -T 0'
        subprocess.run(cmd, shell = True)

#%% Correcting wrong virus names

Name_re = re.compile('\((.+)\)$')

full_Name_re = re.compile('^[AB]/.+/.+/\d{4}.*$')

no_Strain_re = re.compile('^[AB]/([^/]+)/\d{2,4}([^/]*)$')

no_full_Year_re = re.compile('^[AB]/([^/]+)/([^/]+)/\d{2}(.*)$')

Description_re = re.compile('^gi\|[0-9]+\|gb\|.*\|Influenza ([AB]/.+/\d+)[ ,]')

metadata_cluster = {}

for Type in ['A', 'B']:
    metadata_cluster[Type] = {}
    for Segment in list('12345678'):
        metadata_cluster[Type][Segment] = {}
        names = []
        dupe_names = []
        with open('./Data/Clusters/'+Type+'_'+Segment+'.fasta', 'r') as file:
            for seq_record in SeqIO.parse(file, 'fasta'):              
                ID = seq_record.id
                metadata_cluster[Type][Segment][ID] = metadata[Type][Segment][ID]
                Virus_name = metadata_cluster[Type][Segment][ID]['Virus_name']
                if Name_re.search(Virus_name):
                    name = Name_re.search(Virus_name).group(1)
                    if full_Name_re.search(name):
                        Name = name
                    elif no_Strain_re.search(name):
                        Place = no_Strain_re.search(name).group(1)
                        Strain = 'unknown'
                        Year = str(metadata_cluster[Type][Segment][ID]['Year'])
                        Subtype = no_Strain_re.search(name).group(2)
                        Name = Type+'/'+Place+'/'+Strain+'/'+Year+Subtype  
                    elif no_full_Year_re.search(name):
                        Place = no_full_Year_re.search(name).group(1)
                        Strain = no_full_Year_re.search(name).group(2)
                        Year = str(metadata_cluster[Type][Segment][ID]['Year'])
                        Subtype = no_full_Year_re.search(name).group(3)
                        Name = Type+'/'+Place+'/'+Strain+'/'+Year+Subtype
                else:
                    if Type == 'A':
                        Country = metadata_cluster[Type][Segment][ID]['Country']
                        Strain = 'unknown'
                        Year = str(metadata_cluster[Type][Segment][ID]['Year'])
                        Subtype = metadata_cluster[Type][Segment][ID]['Subtype']
                        Name = Type+'/'+Country+'/'+Strain+'/'+Year+'('+Subtype+')'                        
                    if Type == 'B':
                        Description = metadata_cluster[Type][Segment][ID]['Description']
                        if Description_re.search(Description):
                            name = Description_re.search(Description).group(1)
                            if full_Name_re.search(name):
                                Name = name
                            elif no_Strain_re.search(name):
                                Place = no_Strain_re.search(name).group(1)
                                Strain = 'unknown'
                                Year = str(metadata_cluster[Type][Segment][ID]['Year'])
                                Subtype = no_Strain_re.search(name).group(2)
                                Name = Type+'/'+Place+'/'+Strain+'/'+Year+Subtype                                        
                            elif no_full_Year_re.search(name):
                                Place = no_full_Year_re.search(name).group(1)
                                Strain = no_full_Year_re.search(name).group(2)
                                Year = str(metadata_cluster[Type][Segment][ID]['Year'])
                                Name = Type+'/'+Place+'/'+Strain+'/'+Year
                        else:
                            Country = metadata_cluster[Type][Segment][ID]['Country']
                            Strain = 'unknown'
                            Year = str(metadata_cluster[Type][Segment][ID]['Year'])  
                            Name = Type+'/'+Country+'/'+Strain+'/'+Year
                if "'" in Name:
                    Name = Name.replace("'", '"')
                metadata_cluster[Type][Segment][ID]['Name'] = Name
                if Name not in names:
                    names.append(Name)
                else:
                    dupe_names.append(Name)
            for ID in metadata_cluster[Type][Segment]:
                if metadata_cluster[Type][Segment][ID]['Name'] in dupe_names:
                    metadata_cluster[Type][Segment][ID]['Name'] = metadata_cluster[Type][Segment][ID]['Name'] + '_' + ID

# writing the dictionary into .JSON files
    
with open('./Data/Metadata/metadata_cluster.json', 'w') as file:
    json.dump(metadata_cluster, file, default = str)
    
    
#%% Making new fasta file with virus name

for Type in ['A', 'B']:
    for Segment in list('12345678'):
        with open('./Data/Tree/'+Type+'_'+Segment+'.fasta', 'w') as new_fasta:
            with open('./Data/Aligned/'+Type+'_'+Segment+'.fasta', 'r') as f:
                for seq_record in SeqIO.parse(f, 'fasta'):
                    ID = seq_record.id
                    if ID in metadata_cluster[Type][Segment]:
                        seq = str(seq_record.seq)
                        name = metadata_cluster[Type][Segment][ID]['Name']
                        new_fasta.write('>' + name + '\n' + seq + '\n')
 
#%% Creating a tree via FastTree


for Type in ['A', 'B']:
    for Segment in list('12345678'):
        input_file = './Data/Tree/'+Type+'_'+Segment+'.fasta'
        output_file = './Data/Tree/'+Type+'_'+Segment+'.tree'
        cmd = 'fasttree -nt -quote '+input_file+' > '+output_file
        subprocess.run(cmd, shell = True)