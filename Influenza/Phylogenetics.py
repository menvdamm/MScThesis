#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 15 13:32:43 2021

@author: mennovandamme
"""

#%% Dependencies

# conda install -c bioconda cd-hit
# conda install -c bioconda fasttree
import time, os, subprocess, json, re
from Bio import SeqIO

#%% Directories

if not os.path.isdir('./Data/Tree'): os.mkdir('./Data/Tree')
if not os.path.isdir('./Data/Tree/Clustered'): os.mkdir('./Data/Tree/Clustered')

#%% Loading the necessary data

with open('./Data/Metadata/metadata_human.json') as file:
    metadata_human = json.load(file)
    
#%% All clusterings

# with open('Cluster_INFO.txt', 'w') as info_file:
#     for Type in ['A', 'B']:
#         for Segment in list('12345678'):
#             info_file.write('Clustering '+Type+'_'+Segment+':\n')
#             for identity in ['0.995', '0.99', '0.98', '0.97', '0.96', '0.95']:
#                 input_file = './Data/Clean/Human/'+Type+'_'+Segment+'.fasta'
#                 output_file = './Data/Tree/Clustered/'+Type+'_'+Segment+'.fasta'
#                 cmd = 'cd-hit-est -i '+input_file+' -o '+output_file+' -c '+identity+' -n 10 -T 0'
#                 subprocess.run(cmd, shell = True)
#                 cluster_count = 0
#                 seq_counts = []
#                 with open('./Data/Tree/Clustered/'+Type+'_'+Segment+'.fasta.clstr', 'r') as file:
#                     line_list = file.readlines()
#                     for count in range(0, len(line_list)):
#                         line = line_list[count]
#                         if line.startswith('>Cluster'):
#                             cluster_count += 1
#                             if cluster_count >= 1:
#                                 seq_counts.append(line_list[count-1].split()[0])
#                 seq_counts.append(line_list[-1].split()[0])
#                 seq_counts = list(map(int, seq_counts))
#                 info_file.write('\tIdentity: {} \n\t\tclusters: {}\n\t\tmax_seq: {}\n\t\tmin_seq: {}\n'.format(identity, cluster_count, max(seq_counts), min(seq_counts)))

#%% Clustering
t0 = time.time()

for Type in ['A']:
    for Segment in list('12345678'):
        identity = '0.99'
        input_file = './Data/Clean/Human/'+Type+'_'+Segment+'.fasta'
        output_file = './Data/Tree/Clustered/'+Type+'_'+Segment+'_TEMP.fasta'
        cmd = 'cd-hit-est -i '+input_file+' -o '+output_file+' -c '+identity+' -n 10'
        subprocess.run(cmd, shell = True)
        
for Type in ['B']:
    for Segment in list('12345678'):
        identity = '0.995'
        input_file = './Data/Clean/Human/'+Type+'_'+Segment+'.fasta'
        output_file = './Data/Tree/Clustered/'+Type+'_'+Segment+'_TEMP.fasta'
        cmd = 'cd-hit-est -i '+input_file+' -o '+output_file+' -c '+identity+' -n 10'
        subprocess.run(cmd, shell = True)

# renaming the files
for file in os.listdir('./Data/Tree/Clustered/'):
    if file.endswith('.fasta.clstr') and not file.startswith('._'):
        os.rename('./Data/Tree/Clustered/'+file, './Data/Tree/Clustered/'+file[0:3]+'.clstr')

t1 = time.time()
print('Clustering took', round((t1-t0)/60, 1), 'minutes')

#%% Correcting wrong virus names

# finding entries with wrong virus names
#
# Name_re = re.compile('\((.+)\)$')
#
# metadata_cluster = {}
# wrong_names = []
#
# for Type in ['A', 'B']:
#     metadata_cluster[Type] = {}
#     for Segment in list('12345678'):
#         metadata_cluster[Type][Segment] = {}
#         with open('./Data/Tree/Clustered/'+Type+'_'+Segment+'_TEMP.fasta', 'r') as file:
#             records = SeqIO.parse(file, 'fasta')
#             for seq_record in records:
#                 ID = seq_record.id
#                 metadata_cluster[Type][Segment][ID] = metadata_human[Type][Segment][ID]
#                 Virus_name = metadata_human[Type][Segment][ID]['Virus_name']
#                 if not Name_re.search(Virus_name):
#                     wrong_names.append(ID)
#                
# with open('wrong_names.txt', 'w') as file:
#     for ID in wrong_names:
#         file.write(ID+'\n')
#         for key in ['Type', 'Subtype', 'Country', 'Year', 'Virus_name', 'Description']:
#             file.write('\t'+key+': '+str(metadata_clean[ID][key])+'\n')

# manually fixing the virus names of the wrong entries
#
# with open('./Data/Metadata/metadata_clean.json') as file:
#    metadata_clean = json.load(file)
#
# Description_re = re.compile('^gi\|[0-9]+\|gb\|.*\|Influenza B/(.*)/(.*)/\d+')
#
# for ID in wrong_names:
#     Type = metadata_clean[ID]['Type']
#     Segment = metadata_clean[ID]['Segment']
#     if Type == 'A':
#         Country = metadata_clean[ID]['Country']
#         Strain = 'unknown'
#         Year = str(metadata_clean[ID]['Year'])
#         Subtype = metadata_clean[ID]['Subtype']
#         Name = Type+'/'+Country+'/'+Strain+'/'+Year+'('+Subtype+')'
#         metadata_cluster[Type][Segment][ID]['Virus_name'] = 'Influenza A virus ('+Name+')'
#     if Type == 'B':
#         Description = metadata_clean[ID]['Description']
#         if Description_re.search(Description):
#             Place = Description_re.search(Description).group(1)
#             Strain = Description_re.search(Description).group(2)
#             Year = str(metadata_clean[ID]['Year'])
#             Name = Type+'/'+Place+'/'+Strain+'/'+Year
#             metadata_cluster[Type][Segment][ID]['Virus_name'] = 'Influenza B virus ('+Name+')'
#         else:
#             Country = metadata_clean[ID]['Country']
#             Strain = 'unknown'
#             Year = str(metadata_clean[ID]['Year'])  
#             Name = Type+'/'+Country+'/'+Strain+'/'+Year
#             metadata_cluster[Type][Segment][ID]['Virus_name'] = 'Influenza B virus ('+Name+')'

# making a seperate metadata dictionary for the representative cluster sequences
# finding and fixing wrong virus names:
#       missing sample collection locations are replaced by collection country
#       missing strain types are replaced by 'unknown'
#       collection years that are represented by only the last 2 digits (eg 97) are replace by the full year (eg 1997)
#       no accents (') allowed in FastTree algorithm => replaced by double accent (")  
# concatinating duplicate names (per type & segment) with the seq ID
# (duplicate names aren't accepted by the FastTree algorithm)

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
        with open('./Data/Tree/Clustered/'+Type+'_'+Segment+'_TEMP.fasta', 'r') as file:
            records = SeqIO.parse(file, 'fasta')
            for seq_record in records:              
                ID = seq_record.id
                metadata_cluster[Type][Segment][ID] = metadata_human[Type][Segment][ID]
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

## other option: just concatinate all names with their ID

# writing the dictionary into .JSON files
    
with open('./Data/Metadata/metadata_cluster.json', 'w') as file:
    json.dump(metadata_cluster, file, default = str)
    
    
#%% Making new fasta file with virus name

for Type in ['A', 'B']:
    for Segment in list('12345678'):
        new_fasta = './Data/Tree/'+Type+'_'+Segment+'.fasta'
        with open(new_fasta, 'w') as new_fasta:
            for ID in metadata_cluster[Type][Segment]:
                Name = metadata_cluster[Type][Segment][ID]['Name']
                Seq = metadata_cluster[Type][Segment][ID]['Aligned_sequence']
                new_fasta.write('>' + Name + '\n' + Seq + '\n')
        
#%% Creating a tree via FastTree

t0 = time.time()

for Type in ['A', 'B']:
    for Segment in list('12345678'):
        input_file = './Data/Tree/'+Type+'_'+Segment+'.fasta'
        output_file = './Data/Tree/'+Type+'_'+Segment+'.tree'
        cmd = 'fasttree -nt -quote '+input_file+' > '+output_file
        subprocess.run(cmd, shell = True)

t1 = time.time()
print('Tree building took', round((t1-t0)/60, 1), 'minutes')
