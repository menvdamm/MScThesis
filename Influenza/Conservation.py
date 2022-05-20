#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 15 17:54:48 2021

@author: mennovandamme
"""

### Creating a dataframe for all data

#%% Dependencies

import json, os
import pandas as pd
import numpy as np
from Bio import SeqIO

#%% Directories

if not os.path.isdir('./Data/Dataframes'): os.mkdir('./Data/Dataframes')
if not os.path.isdir('./Data/Dataframes/Complete'): os.mkdir('./Data/Dataframes/Complete')
if not os.path.isdir('./Data/Dataframes/Ungapped'): os.mkdir('./Data/Dataframes/Ungapped')

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
        tot_length += metadata[ID]['Seq_length']    
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
    print('\t \t Complete df length: '+len(complete_df)+'\n',
          '\t \t Average length: '+avg_length+'\n'
          '\t \t Ungapped df length: '+df_length+'\n',
          '\t \t Allowed Gap Percentage: '+allowed_gap_percentage)
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

#%% Dataframe per type





#%% Complete dataframe per segment

# def complete_segment_df_maker(Type, Segment):
#     oldest_seq = min(metadata_human[Type][Segment].items(), key = lambda item: item[1]['Year'])[1]['Aligned_sequence']
#     seq_length = len(oldest_seq)
#     Positions = list(range(1, seq_length + 1))
#     Mutation_count = [0] * seq_length
#     Variability = [0] * seq_length
#     Entropy = [0] * seq_length
#     total_seq_count = 0
#     Gap_count = [0] * seq_length
#     A_count = [0] * seq_length
#     T_count = [0] * seq_length
#     G_count = [0] * seq_length
#     C_count = [0] * seq_length
#     Date_ID = {}
#     for ID in metadata_human[Type][Segment]:
#         Year = metadata_human[Type][Segment][ID]['Year']
#         Month =  metadata_human[Type][Segment][ID]['Month']
#         Day =  metadata_human[Type][Segment][ID]['Day']
#         if not (Year == 'unknown' or Month == 'unknown' or Day == 'unknown'):
#             if Year not in Date_ID.keys():
#                 Date_ID[Year] = {}
#             if Month not in Date_ID[Year].keys():
#                 Date_ID[Year][Month] = {}
#             if Day not in Date_ID[Year][Month].keys():
#                 Date_ID[Year][Month][Day] = [ID]
#             else:
#                 Date_ID[Year][Month][Day].append(ID)
#     previous_seq = oldest_seq
#     for Year in sorted(Date_ID.keys()):
#         yearly_seq_count = sum([len(Date_ID[Year][Month][Day]) for Month in Date_ID[Year] for Day in Date_ID[Year][Month]])          
#         for Month in sorted(Date_ID[Year].keys()):
#             for Day in sorted(Date_ID[Year][Month].keys()):
#                 for ID in Date_ID[Year][Month][Day]:
#                     total_seq_count += 1
#                     seq = metadata_human[Type][Segment][ID]['Aligned_sequence']
#                     for pos in Positions:
#                         if seq[pos-1] != previous_seq[pos-1]:
#                             Mutation_count[pos-1] += 1
#                             Variability[pos-1] += 1/yearly_seq_count
#                         if seq[pos-1] == '-':
#                             Gap_count[pos-1] += 1
#                         elif seq[pos-1] == 'A':
#                             A_count[pos-1] += 1
#                         elif seq[pos-1] == 'T' or seq[pos-1] == 'U':
#                             T_count[pos-1] += 1
#                         elif seq[pos-1] == 'G':
#                             G_count[pos-1] += 1
#                         elif seq[pos-1] == 'C':
#                             C_count[pos-1] += 1  
#                     previous_seq = seq
#     Gap_ratio = np.divide(Gap_count, total_seq_count)
#     A_ratio = np.divide(A_count, total_seq_count)
#     T_ratio = np.divide(T_count, total_seq_count)
#     G_ratio = np.divide(G_count, total_seq_count)
#     C_ratio = np.divide(C_count, total_seq_count)
#     for i in range(0, seq_length):
#         if Gap_ratio[i] == 0:
#             E_gap = 0
#         else:
#             E_gap = Gap_ratio[i] * np.log2(Gap_ratio[i])
#         if A_ratio[i] == 0:
#             E_A = 0        
#         else:
#             E_A = A_ratio[i] * np.log2(A_ratio[i])
#         if T_ratio[i] == 0:
#             E_T = 0
#         else:
#             E_T = T_ratio[i] * np.log2(T_ratio[i])
#         if G_ratio[i] == 0:
#             E_G = 0 
#         else:
#             E_G = G_ratio[i] * np.log2(G_ratio[i])
#         if C_ratio[i] == 0:
#             E_C = 0 
#         else: 
#             E_C = C_ratio[i] * np.log2(C_ratio[i])
#         Entropy[i] = abs(- E_gap - E_A - E_T - E_G - E_C)
#     df = pd.DataFrame({'Position': Positions,
#                        'Mutation_count': Mutation_count,
#                        'Variability': Variability,
#                        'Shannon_entropy': Entropy,
#                        'Gap_percentage': Gap_ratio * 100,
#                        'A_percentage': A_ratio * 100,
#                        'T_percentage': T_ratio * 100,
#                        'G_percentage': G_ratio * 100,
#                        'C_percentage': C_ratio * 100})
#     df['Position'] = df['Position'].astype('int64')
#     df['Mutation_count'] = df['Mutation_count'].astype('int64')
#     return df

#%% Complete dataframe per type

# for Type in ['A','B']:
#     for Segment in list('12345678'):
#         with open('./Data/Dataframes/Complete/complete_df_'+Type+'_'+Segment+'.csv') as file:
#             globals()['complete_df_'+Type+'_'+Segment] = pd.read_csv(file)


# def complete_type_df_maker(Type):
#     df = pd.DataFrame({'Position': [],
#                        'Segment': [],
#                        'Segment_position': [],
#                        'Mutation_count': [],
#                        'Variability': [],
#                        'Shannon_entropy': [],
#                        'Gap_percentage': [],
#                        'A_percentage': [],
#                        'T_percentage': [],
#                        'G_percentage': [],
#                        'C_percentage': []})
#     length = 0
#     for Segment in list('12345678'):
#         df_Segment = globals()['complete_df_'+Type+'_'+Segment].copy()
#         df_Segment['Segment'] = [int(Segment)] * len(df_Segment)
#         df_Segment = df_Segment.rename({'Position': 'Segment_position'}, axis='columns')
#         df_Segment['Position']  = df_Segment['Segment_position'] + length
#         length = max(df_Segment['Position'])
#         df = df.append(df_Segment)
#     df['Position'] = df['Position'].astype('int64')
#     df['Segment'] = df['Segment'].astype('int64')
#     df['Segment_position'] = df['Segment_position'].astype('int64')
#     df['Mutation_count'] = df['Mutation_count'].astype('int64')
#     return df

# for Type in ['A','B']:
#     globals()['complete_df_'+Type] = complete_type_df_maker(Type)
#     globals()['complete_df_'+Type].to_csv('./Data/Dataframes/Complete/complete_df_'+Type+'.csv', index = False)   

#%% Ungapped dataframe per segment

# def segment_df_ungapper(Type, Segment):
#     df = globals()['complete_df_'+Type+'_'+Segment].copy()
#     df = df[df['Gap_percentage'] < 90].sort_values('Position')
#     df = df.rename({'Position': 'Aligned_position'}, axis='columns')
#     df['Position'] = list(range(1, len(df)+1))
#     return df

# for Type in ['A', 'B']:
#     for Segment in list('12345678'):
#         globals()['df_'+Type+'_'+Segment] = segment_df_ungapper(Type, Segment)
#         globals()['df_'+Type+'_'+Segment].to_csv('./Data/Dataframes/Ungapped/df_'+Type+'_'+Segment+'.csv', index = False)
        
#%% Ungapped dataframe per type

# def type_df_ungapper(Type):
#     df = globals()['complete_df_'+Type].copy()
#     df = df[df['Gap_percentage'] < 90].sort_values('Position')
#     df = df.rename({'Position': 'Aligned_position', 'Segment_position': 'Aligned_segment_position'}, axis='columns')
#     df['Position'] = list(range(1, len(df)+1))
#     Segment_position = []
#     for Segment in list('12345678'):
#         df_Segment = df[df['Segment'] == int(Segment)]
#         Segment_position += list(range(1, len(df_Segment)+1))
#     df['Segment_position'] = Segment_position
#     return df

# for Type in ['A', 'B']:
#     globals()['df_'+Type] = type_df_ungapper(Type)
#     globals()['df_'+Type].to_csv('./Data/Dataframes/Ungapped/df_'+Type+'.csv', index = False)
        