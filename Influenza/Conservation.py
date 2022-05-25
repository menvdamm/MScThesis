#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 15 17:54:48 2021

@author: mennovandamme
"""

### Creating a dataframe for all data

#%% Dependencies

import json, os, subprocess, re
import pandas as pd
import numpy as np
from Bio import SeqIO

#%% Directories

if not os.path.isdir('./Data/Dataframes'): os.mkdir('./Data/Dataframes')
if not os.path.isdir('./Data/Dataframes/Complete'): os.mkdir('./Data/Dataframes/Complete')
if not os.path.isdir('./Data/Dataframes/Ungapped'): os.mkdir('./Data/Dataframes/Ungapped')
if not os.path.isdir('./Data/Dataframes/Scores'): os.mkdir('./Data/Dataframes/Scores')
if not os.path.isdir('./Data/Genome'): os.mkdir('./Data/Genome')
if not os.path.isdir('./Data/CDS'): os.mkdir('./Data/CDS')
if not os.path.isdir('./Data/Consensus'): os.mkdir('./Data/Consensus')

#%% Loading the necessary data

with open('./Data/Metadata/metadata.json') as file:
    metadata = json.load(file)

#%% Function to sort IDs by date

def sort_IDs_by_date(IDs, metadata):
    date_IDs = {}
    IDs_sorted = []    
    if len(IDs) == 1:
        return IDs
    else:
        for ID in IDs:
            Year = metadata[ID]['Year']
            Month = metadata[ID]['Month']
            Day = metadata[ID]['Day']
            if Year not in date_IDs:
                date_IDs[Year] = {}
            if Month not in date_IDs[Year]:
                date_IDs[Year][Month] = {}
            if Day not in date_IDs[Year][Month]:
                date_IDs[Year][Month][Day] = []
            date_IDs[Year][Month][Day].append(ID)
        for Year in sorted(date_IDs.keys()):
            for Month in sorted(date_IDs[Year].keys()):
                for Day in sorted(date_IDs[Year][Month].keys()):
                    for ID in date_IDs[Year][Month][Day]:
                        IDs_sorted.append(ID)  
    return IDs_sorted

#%% Complete dataframe per segment

def make_df(fasta, metadata):
    seq_dict = SeqIO.index(fasta, 'fasta')
    IDs_sorted = sort_IDs_by_date(list(metadata.keys()), metadata)
    oldest_seq = str(seq_dict[IDs_sorted[0]].seq)
    seq_length = len(oldest_seq)
    Positions = list(range(1, seq_length + 1))
    seq_count = 0
    Mut_count = [0] * seq_length
    Entropy = [0] * seq_length
    Gap_count = [0] * seq_length
    A_count = [0] * seq_length
    T_count = [0] * seq_length
    G_count = [0] * seq_length
    C_count = [0] * seq_length
    previous_seq = oldest_seq
   
    for ID in IDs_sorted:
        seq_record = seq_dict[ID]
        seq = str(seq_record.seq)
        seq_count += 1
        for i in range(seq_length):
            if seq[i] != previous_seq[i]:
                Mut_count[i] += 1
            if seq[i] == '-':
                Gap_count[i] += 1
            elif seq[i] == 'A':
                A_count[i] += 1
            elif seq[i] == 'T':
                T_count[i] += 1
            elif seq[i] == 'G':
                G_count[i] += 1
            elif seq[i] == 'C':
                C_count[i] += 1  
        previous_seq = seq
        
    Mutability = np.divide(Mut_count, seq_count)
    Gap_ratio = np.divide(Gap_count, seq_count)
    A_ratio = np.divide(A_count, seq_count)
    T_ratio = np.divide(T_count, seq_count)
    G_ratio = np.divide(G_count, seq_count)
    C_ratio = np.divide(C_count, seq_count)
    
    for i in range(seq_length):       
        if Gap_ratio[i] == 0:
            E_gap = 0
        else:
            E_gap = Gap_ratio[i] * np.log2(Gap_ratio[i])
        if A_ratio[i] == 0:
            E_A = 0        
        else:
            E_A = A_ratio[i] * np.log2(A_ratio[i])
        if T_ratio[i] == 0:
            E_T = 0
        else:
            E_T = T_ratio[i] * np.log2(T_ratio[i])
        if G_ratio[i] == 0:
            E_G = 0 
        else:
            E_G = G_ratio[i] * np.log2(G_ratio[i])
        if C_ratio[i] == 0:
            E_C = 0 
        else: 
            E_C = C_ratio[i] * np.log2(C_ratio[i])           
        Entropy[i] = abs(- E_gap - E_A - E_T - E_G - E_C)
        
    df = pd.DataFrame({'Position': Positions,
                       'Mutability': Mutability,
                       'Shannon_entropy': Entropy,
                       'Gap_percentage': Gap_ratio * 100,
                       'A_percentage': A_ratio * 100,
                       'T_percentage': T_ratio * 100,
                       'G_percentage': G_ratio * 100,
                       'C_percentage': C_ratio * 100})    
    df['Position'] = df['Position'].astype('int64')
    return df

for Type in ['A','B']:
    for Segment in list('12345678'):
        fasta_file = './Data/Aligned/'+Type+'_'+Segment+'.fasta'
        globals()['complete_df_'+Type+'_'+Segment] = make_df(fasta_file, metadata[Type][Segment])
        globals()['complete_df_'+Type+'_'+Segment].to_csv('./Data/Dataframes/Complete/df_'+Type+'_'+Segment+'.csv', index = False)   

#%% Ungapped dataframe per segment

# Identifying rare insertions, allowed gap percentage
def find_rare_insertions(complete_df, metadata):
    tot_length = 0
    for ID in metadata:
        tot_length += int(metadata[ID]['Seq_length'])
    avg_length = int(tot_length/len(metadata)) 
    for allowed_gap_percentage in list(range(100, 0, -1)):
        smaller = complete_df['Gap_percentage'] < allowed_gap_percentage
        df_length = smaller.sum()
        if df_length <= avg_length:
            break  
    return avg_length, df_length, allowed_gap_percentage

# removing the positions from the complete df that have a gap percentage over the allowed gap percentage
def ungap_df(complete_df, metadata):
    avg_length, df_length, allowed_gap_percentage = find_rare_insertions(complete_df, metadata)
    print('\t\t Complete df length: '+str(len(complete_df.index))+'\n',
          '\t\t Average length: '+str(avg_length)+'\n'
          '\t\t Ungapped df length: '+str(df_length)+'\n',
          '\t\t Allowed Gap Percentage: '+str(allowed_gap_percentage))
    df = complete_df
    df = df[df['Gap_percentage'] < allowed_gap_percentage].sort_values('Position')
    df['Position'] = list(range(1, len(df)+1))
    return df

for Type in ['A','B']:
    print('Type: '+Type)
    for Segment in list('12345678'):
        print('\t Segment: '+Segment)
        complete_df = globals()['complete_df_'+Type+'_'+Segment] 
        globals()['df_'+Type+'_'+Segment] = ungap_df(complete_df, metadata[Type][Segment])
        globals()['df_'+Type+'_'+Segment].to_csv('./Data/Dataframes/Ungapped/df_'+Type+'_'+Segment+'.csv', index = False)

# Type: A
# 	 Segment: 1
# 		 Complete df length: 2692
#  		 Average length: 2303
# 		 Ungapped df length: 2297
#  		 Allowed Gap Percentage: 38
# 	 Segment: 2
# 		 Complete df length: 2447
#  		 Average length: 2303
# 		 Ungapped df length: 2303
#  		 Allowed Gap Percentage: 34
# 	 Segment: 3
# 		 Complete df length: 2306
#  		 Average length: 2189
# 		 Ungapped df length: 2189
#  		 Allowed Gap Percentage: 39
# 	 Segment: 4
# 		 Complete df length: 2283
#  		 Average length: 1724
# 		 Ungapped df length: 1718
#  		 Allowed Gap Percentage: 52
# 	 Segment: 5
# 		 Complete df length: 1817
#  		 Average length: 1527
# 		 Ungapped df length: 1526
#  		 Allowed Gap Percentage: 38
# 	 Segment: 6
# 		 Complete df length: 1877
#  		 Average length: 1426
# 		 Ungapped df length: 1426
#  		 Allowed Gap Percentage: 51
# 	 Segment: 7
# 		 Complete df length: 1266
#  		 Average length: 994
# 		 Ungapped df length: 994
#  		 Allowed Gap Percentage: 48
# 	 Segment: 8
# 		 Complete df length: 990
#  		 Average length: 857
# 		 Ungapped df length: 856
#  		 Allowed Gap Percentage: 40
# Type: B
# 	 Segment: 1
# 		 Complete df length: 2796
#  		 Average length: 2317
# 		 Ungapped df length: 2312
#  		 Allowed Gap Percentage: 27
# 	 Segment: 2
# 		 Complete df length: 2802
#  		 Average length: 2349
# 		 Ungapped df length: 2349
#  		 Allowed Gap Percentage: 29
# 	 Segment: 3
# 		 Complete df length: 2317
#  		 Average length: 2250
# 		 Ungapped df length: 2250
#  		 Allowed Gap Percentage: 27
# 	 Segment: 4
# 		 Complete df length: 1909
#  		 Average length: 1812
# 		 Ungapped df length: 1798
#  		 Allowed Gap Percentage: 40
# 	 Segment: 5
# 		 Complete df length: 1857
#  		 Average length: 1784
# 		 Ungapped df length: 1776
#  		 Allowed Gap Percentage: 23
# 	 Segment: 6
# 		 Complete df length: 1744
#  		 Average length: 1501
# 		 Ungapped df length: 1495
#  		 Allowed Gap Percentage: 22
# 	 Segment: 7
# 		 Complete df length: 1197
#  		 Average length: 1142
# 		 Ungapped df length: 1142
#  		 Allowed Gap Percentage: 26
# 	 Segment: 8
# 		 Complete df length: 1122
#  		 Average length: 1051
# 		 Ungapped df length: 1050
#  		 Allowed Gap Percentage: 39
          
#%% Dataframe per type

def make_type_df(Type):
    df = pd.DataFrame({'Position': [],
                        'Segment': [],
                        'Segment_position': [],
                        'Mutability': [],
                        'Shannon_entropy': [],
                        'Gap_percentage': [],
                        'A_percentage': [],
                        'T_percentage': [],
                        'G_percentage': [],
                        'C_percentage': []})
    length = 0
    for Segment in list('12345678'):
        df_Segment = globals()['df_'+Type+'_'+Segment].copy()
        df_Segment['Segment'] = [int(Segment)] * len(df_Segment)
        df_Segment = df_Segment.rename({'Position': 'Segment_position'}, axis='columns')
        df_Segment['Position']  = df_Segment['Segment_position'] + length
        length = max(df_Segment['Position'])
        df = df.append(df_Segment)
    df['Position'] = df['Position'].astype('int64')
    df['Segment'] = df['Segment'].astype('int64')
    df['Segment_position'] = df['Segment_position'].astype('int64')
    return df.sort_values('Position')

for Type in ['A','B']:
    globals()['df_'+Type] = make_type_df(Type)
    globals()['df_'+Type].to_csv('./Data/Dataframes/df_'+Type+'.csv', index = False)   

#%% Scoring 30 nucleotide regions per segment

def score_window(df, size):
    Begin = pd.Series(list(range(0, len(df)-size+1)))
    End = pd.Series(list(range(size-1, len(df))))
    Mut_scores = pd.Series([sum(df['Mutability'][begin:end+1]) for begin, end in zip(Begin, End)])
    Ent_scores = pd.Series([sum(df['Shannon_entropy'][begin:end+1]) for begin, end in zip(Begin, End)])
    Gap_scores = pd.Series([sum(df['Gap_percentage'][begin:end+1]) for begin, end in zip(Begin, End)])
    score_df = pd.DataFrame({'Begin_position': Begin + 1,
                             'End_position': End + 1,     
                             'Mutability': Mut_scores,
                             'Shannon_entropy': Ent_scores,
                             'Gap_percentage': Gap_scores})
    score_df = score_df.sort_values('Mutability')
    return score_df  

for Type in ['A','B']:
    for Segment in list('12345678'):
        globals()['score_df_'+Type+'_'+Segment] = score_window(globals()['df_'+Type+'_'+Segment], 30)
        globals()['score_df_'+Type+'_'+Segment].to_csv('./Data/Dataframes/Scores/df_'+Type+'_'+Segment+'.csv', index = False)
        
#%% Score dataframe per type        
       
# makes sure that there are no 30bp regions scored over 2 segments

def make_type_score_df(Type):
    df = pd.DataFrame({'Begin_position': [],
                       'End_position': [],
                        'Segment': [],
                        'Segment_begin_position': [],
                        'Segment_end_position': [],
                        'Mutability': [],
                        'Shannon_entropy': [],
                        'Gap_percentage': []})
    length = 0
    for Segment in list('12345678'):
        df_Segment = globals()['score_df_'+Type+'_'+Segment].copy()
        df_Segment['Segment'] = [int(Segment)] * len(df_Segment)
        df_Segment = df_Segment.rename({'Begin_position': 'Segment_begin_position'}, axis='columns')
        df_Segment = df_Segment.rename({'End_position': 'Segment_end_position'}, axis='columns')
        df_Segment['Begin_position']  = df_Segment['Segment_begin_position'] + length
        df_Segment['End_position']  = df_Segment['Segment_end_position'] + length
        length = max(df_Segment['End_position'])
        df = df.append(df_Segment)
    df['Begin_position'] = df['Begin_position'].astype('int64')
    df['Begin_position'] = df['Begin_position'].astype('int64')
    df['Segment'] = df['Segment'].astype('int64')
    df['Segment_begin_position'] = df['Segment_begin_position'].astype('int64')
    df['Segment_end_position'] = df['Segment_end_position'].astype('int64')
    return df.sort_values('Mutability')

for Type in ['A','B']:
    globals()['complete_score_df_'+Type] = make_type_score_df(Type)
    globals()['complete_score_df_'+Type].to_csv('./Data/Dataframes/complete_score_df_'+Type+'.csv', index = False)  

#%% Making consensus sequence 

# for Type in ['A','B']:
#     with open('./Data/Genome/'+Type+'_consensus.fasta', 'w') as file:
#         consensus = ''
#         for row in globals()['df_'+Type].sort_values('Position').iterrows():
#             percentages = {'A': row[1]['A_percentage'],
#                            'T': row[1]['T_percentage'],
#                            'G': row[1]['G_percentage'],
#                            'C': row[1]['C_percentage']}
#             most_occuring = max(percentages, key = percentages.get)
#             consensus += most_occuring
#         file.write('>' + Type + '\n' + consensus)
        
for Type in ['A','B']:
    for Segment in list('12345678'):
        with open('./Data/Genome/Consensus/'+Type+'_'+Segment+'.fasta', 'w') as file:
            consensus = ''
            for row in globals()['df_'+Type+'_'+Segment].sort_values('Position').iterrows():
                percentages = {'A': row[1]['A_percentage'],
                               'T': row[1]['T_percentage'],
                               'G': row[1]['G_percentage'],
                               'C': row[1]['C_percentage']}
                most_occuring = max(percentages, key = percentages.get)
                consensus += most_occuring
            file.write('>' + Type + '_' + Segment + '\n' + consensus)  
        
#%% Align CDS to consensus

# manually downloaded coding sequences from NCBI Influenza virus reference sequence
# https://www.ncbi.nlm.nih.gov/nuccore/1798174254?report=genbank
# => './Data/Genome/A_CDS.fasta' (M2 & NEP gene had to be split for proper alignment)
# => './Data/Genome/B_CDS.fasta' (NEP gene had to be split for proper alignment)

for Type in ['A','B']:
     for Segment in list('12345678'):
         CDS_file = './Data/Genome/'+Type+'_CDS.fasta'
         consensus_file = './Data/Genome/'+Type+'_consensus.fasta'
         output_file = './Data/Genome/'+Type+'.fasta'
         cmd = 'mafft --auto  --preservecase --keeplength --addfragments '+CDS_file+' --reorder --thread -1 '+consensus_file+' > '+output_file
         subprocess.run(cmd, shell = True) 


# for Type in ['A','B']:
#     for Segment in list('12345678'):
#         CDS_file = './Data/Genome/'+Type+'_CDS.fasta'
#         consensus_file = './Data/Genome/'+Type+'_consensus.fasta'
#         output_file = './Data/Genome/'+Type+'.fasta'
#         cmd = 'mafft --auto  --preservecase --keeplength --addfragments '+CDS_file+' --reorder --thread -1 '+consensus_file+' > '+output_file
#         subprocess.run(cmd, shell = True)   

# CDS_file = './Data/Genome/A_HA.fasta'
# consensus_file = './Data/Genome/A_4_consensus.fasta'
# output_file = './Data/Genome/A_4_aligned.fasta'
# cmd = 'mafft --auto  --preservecase --keeplength --add '+CDS_file+' --reorder --thread -1 '+consensus_file+' > '+output_file
# subprocess.run(cmd, shell = True)   

# input_file = './Data/Genome/A_4_combined.fasta'
# output_file = './Data/Genome/A_4_aligned.fasta'
# cmd = 'mafft --thread -1 --auto --preservecase '+input_file+' > '+output_file
# subprocess.run(cmd, shell = True)    

#%% Extract positions of coding sequences

protein_re = re.compile('\[gene=([^\]]*)\]')

for Type in ['A','B']:
    begin_pos = []
    end_pos = []
    proteins = []
    for seq_record in SeqIO.parse('./Data/Genome/'+Type+'.fasta', 'fasta'):
        ID = seq_record.id
        if ID != Type:
            protein = protein_re.search(seq_record.description).group(1)
            proteins.append(protein)
            seq = str(seq_record.seq)
            for i in range(len(seq)):
                 if seq[i] != '-':
                     begin_pos.append(i+1)
                     break
            for i in range(len(seq)-1,0,-1):
                 if seq[i] != '-':
                     end_pos.append(i+1)
                     break            
    globals()['CDS_df_'+Type] = pd.DataFrame({'Begin_position': begin_pos,
                                               'End_position': end_pos,
                                               'Protein': proteins})                
    globals()['CDS_df_'+Type].to_csv('./Data/Genome/CDS_df_'+Type+'.csv', index = False)        
        
#%% Filter out any 30bp regions that aren't completely in the coding regions

def filter_scores(df, CDS_df):
    coding_regions = []
    for i in range(len(CDS_df)):
        b = CDS_df['Begin_position'][i]
        e = CDS_df['End_position'][i]
        coding_regions += list(range(b, e + 1))
    coding_regions = set(coding_regions)
    df = df[df["Begin_position"].isin(coding_regions)]
    score_df = df[df["End_position"].isin(coding_regions)]
    return score_df  

for Type in ['A','B']:
    globals()['score_df_'+Type] = filter_scores(globals()['complete_score_df_'+Type], globals()['CDS_df_'+Type])
    globals()['score_df_'+Type].to_csv('./Data/Dataframes/score_df_'+Type+'.csv', index = False)        
        
        
        
        
        
        
        
        