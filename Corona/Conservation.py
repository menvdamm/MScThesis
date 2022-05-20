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

#%% Loading the necessary data

with open('./Data/Metadata/metadata.json', 'r') as file:
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

complete_df = make_df('./Data/Alignment/SARSCoV2.fasta', metadata)
complete_df.to_csv('./Data/Dataframes/complete_df.csv', index = False)
    
#%% Identifying rare insertions, allowed gap percentage

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

avg_length, df_length, allowed_gap_percentage = find_rare_insertions(complete_df, metadata)
# 29784, 29783, 45

#%% Ungapped dataframe per segment

def ungap_df(df, allowed_gap_percentage):
    df = df[df['Gap_percentage'] < allowed_gap_percentage].sort_values('Position')
    df['Position'] = list(range(1, len(df)+1))
    return df

df = ungap_df(complete_df, allowed_gap_percentage)
df.to_csv('./Data/Dataframes/df.csv', index = False)

#%% Scoring 39 nucleotide regions

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

score_df = score_window(df, 39)
score_df.to_csv('./Data/Dataframes/score_df.csv', index = False)

