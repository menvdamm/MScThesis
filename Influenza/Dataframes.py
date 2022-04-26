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

#%% Directories

if not os.path.isdir('./Data/Dataframes'): os.mkdir('./Data/Dataframes')
if not os.path.isdir('./Data/Dataframes/Animations'): os.mkdir('./Data/Dataframes/Animations')
if not os.path.isdir('./Data/Dataframes/Complete'): os.mkdir('./Data/Dataframes/Complete')
if not os.path.isdir('./Data/Dataframes/Ungapped'): os.mkdir('./Data/Dataframes/Ungapped')
if not os.path.isdir('./Data/Dataframes/Scores'): os.mkdir('./Data/Dataframes/Scores')

#%% Loading the necessary data

with open('./Data/Metadata/metadata_human.json') as file:
    metadata_human = json.load(file)
    
with open('./Data/Metadata/references.json') as file:
    references = json.load(file)

#%% Animation dataframe per segment

def segment_animation_df_maker(Type, Segment):
    oldest_seq = min(metadata_human[Type][Segment].items(), key = lambda item: item[1]['Year'])[1]['Aligned_sequence']
    Positions = list(range(1, len(oldest_seq) + 1))
    Mutation_count = [0] * len(oldest_seq)
    Variability = [0] * len(oldest_seq)
    df = pd.DataFrame({'Position': [],
                       'Mutation_count': [],
                       'Variability': [],
                       'Year': [] })
    Date_ID = {}
    for ID in metadata_human[Type][Segment]:
        Year = metadata_human[Type][Segment][ID]['Year']
        Month =  metadata_human[Type][Segment][ID]['Month']
        Day =  metadata_human[Type][Segment][ID]['Day']
        if not (Year == 'unknown' or Month == 'unknown' or Day == 'unknown'):
            if Year not in Date_ID.keys():
                Date_ID[Year] = {}
            if Month not in Date_ID[Year].keys():
                Date_ID[Year][Month] = {}
            if Day not in Date_ID[Year][Month].keys():
                Date_ID[Year][Month][Day] = [ID]
            else:
                Date_ID[Year][Month][Day].append(ID)
    previous_seq = oldest_seq
    for Year in sorted(Date_ID.keys()):
        Years = [Year] * len(oldest_seq)
        seq_count = sum([len(Date_ID[Year][Month][Day]) for Month in Date_ID[Year] for Day in Date_ID[Year][Month]])          
        for Month in sorted(Date_ID[Year].keys()):
            for Day in sorted(Date_ID[Year][Month].keys()):
                for ID in Date_ID[Year][Month][Day]:        
                    seq = metadata_human[Type][Segment][ID]['Aligned_sequence']
                    for pos in Positions:
                        if seq[pos-1] != previous_seq[pos-1]:
                            Mutation_count[pos-1] += 1
                            Variability[pos-1] += 1/seq_count
                    previous_seq = seq
        df_Year = pd.DataFrame({'Position': Positions,
                                'Mutation_count': Mutation_count,
                                'Variability': Variability,
                                'Year': Years})
        df = df.append(df_Year)
    df['Year'] = df['Year'].astype('int64')
    df['Position'] = df['Position'].astype('int64')
    df['Mutation_count'] = df['Mutation_count'].astype('int64')
    return df

for Type in ['A','B']:
    for Segment in list('12345678'):
        globals()['animation_df_'+Type+'_'+Segment] = segment_animation_df_maker(Type, Segment)
        globals()['animation_df_'+Type+'_'+Segment].to_csv('./Data/Dataframes/Animations/animation_df_'+Type+'_'+Segment+'.csv', index = False)   

#%% Animation dataframe per type

def sub_segment_animation_df_maker(Type, Segment, all_Years):
    oldest_seq = min(metadata_human[Type][Segment].items(), key = lambda item: item[1]['Year'])[1]['Aligned_sequence']
    Positions = list(range(1, len(oldest_seq) + 1))
    Mutation_count = [0] * len(oldest_seq)
    Variability = [0] * len(oldest_seq)
    df = pd.DataFrame({'Position': [],
                       'Mutation_count': [],
                       'Variability': [],
                       'Year': [] })
    Date_ID = {}
    for ID in metadata_human[Type][Segment]:
        Year = metadata_human[Type][Segment][ID]['Year']
        Month =  metadata_human[Type][Segment][ID]['Month']
        Day =  metadata_human[Type][Segment][ID]['Day']
        if not (Year == 'unknown' or Month == 'unknown' or Day == 'unknown'):
            if Year not in Date_ID.keys():
                Date_ID[Year] = {}
            if Month not in Date_ID[Year].keys():
                Date_ID[Year][Month] = {}
            if Day not in Date_ID[Year][Month].keys():
                Date_ID[Year][Month][Day] = [ID]
            else:
                Date_ID[Year][Month][Day].append(ID)
    previous_seq = oldest_seq
    for Year in sorted(all_Years):
        Years = [Year] * len(oldest_seq)
        if Year in Date_ID.keys(): 
            seq_count = sum([len(Date_ID[Year][Month][Day]) for Month in Date_ID[Year] for Day in Date_ID[Year][Month]])          
            for Month in sorted(Date_ID[Year].keys()):
                for Day in sorted(Date_ID[Year][Month].keys()):
                    for ID in Date_ID[Year][Month][Day]:        
                        seq = metadata_human[Type][Segment][ID]['Aligned_sequence']
                        for pos in Positions:
                            if seq[pos-1] != previous_seq[pos-1]:
                                Mutation_count[pos-1] += 1
                                Variability[pos-1] += 1/seq_count
                        previous_seq = seq
            df_Year = pd.DataFrame({'Position': Positions,
                                    'Mutation_count': Mutation_count,
                                    'Variability': Variability,
                                    'Year': Years})
            df = df.append(df_Year)
        elif Year not in Date_ID.keys():
            df_Year = pd.DataFrame({'Position': Positions,
                                    'Mutation_count': Mutation_count,
                                    'Variability': Variability,
                                    'Year': Years})
            df = df.append(df_Year)  
    return df


def type_animation_df_maker(Type):
    df = pd.DataFrame({'Position': [],
                       'Mutation_count': [],
                       'Variability': [],
                       'Year': []})
    all_Years = []
    for Segment in list('12345678'):
        for ID in metadata_human[Type][Segment]:
            Year = metadata_human[Type][Segment][ID]['Year']
            if Year not in all_Years:
                all_Years.append(Year)   
    length = 0
    for Segment in list('12345678'):
        df_Segment = sub_segment_animation_df_maker(Type, Segment, all_Years)
        df_Segment['Position']  = df_Segment['Position'] + length
        length = max(df_Segment['Position'])
        df = df.append(df_Segment)
    df['Year'] = df['Year'].astype('int64')
    df['Position'] = df['Position'].astype('int64')
    df['Mutation_count'] = df['Mutation_count'].astype('int64')
    return df

for Type in ['A','B']:
    globals()['animation_df_'+Type] = type_animation_df_maker(Type)
    globals()['animation_df_'+Type].to_csv('./Data/Dataframes/Animations/animation_df_'+Type+'.csv', index = False)   
    
#%% Complete dataframe per segment

def complete_segment_df_maker(Type, Segment):
    df = pd.DataFrame({'Position': [],
                       'Mutation_count': [],
                       'Variability': [],
                       'Shannon_entropy': [],
                       'Gap_percentage': [],
                       'A_percentage': [],
                       'T_percentage': [],
                       'G_percentage': [],
                       'C_percentage': []})
    oldest_seq = min(metadata_human[Type][Segment].items(), key = lambda item: item[1]['Year'])[1]['Aligned_sequence']
    seq_length = len(oldest_seq)
    Positions = list(range(1, seq_length + 1))
    Mutation_count = [0] * seq_length
    Variability = [0] * seq_length
    Entropy = [0] * seq_length
    total_seq_count = 0
    Gap_count = [0] * seq_length
    A_count = [0] * seq_length
    T_count = [0] * seq_length
    G_count = [0] * seq_length
    C_count = [0] * seq_length
    Date_ID = {}
    for ID in metadata_human[Type][Segment]:
        Year = metadata_human[Type][Segment][ID]['Year']
        Month =  metadata_human[Type][Segment][ID]['Month']
        Day =  metadata_human[Type][Segment][ID]['Day']
        if not (Year == 'unknown' or Month == 'unknown' or Day == 'unknown'):
            if Year not in Date_ID.keys():
                Date_ID[Year] = {}
            if Month not in Date_ID[Year].keys():
                Date_ID[Year][Month] = {}
            if Day not in Date_ID[Year][Month].keys():
                Date_ID[Year][Month][Day] = [ID]
            else:
                Date_ID[Year][Month][Day].append(ID)
    previous_seq = oldest_seq
    for Year in sorted(Date_ID.keys()):
        yearly_seq_count = sum([len(Date_ID[Year][Month][Day]) for Month in Date_ID[Year] for Day in Date_ID[Year][Month]])          
        for Month in sorted(Date_ID[Year].keys()):
            for Day in sorted(Date_ID[Year][Month].keys()):
                for ID in Date_ID[Year][Month][Day]:
                    total_seq_count += 1
                    seq = metadata_human[Type][Segment][ID]['Aligned_sequence']
                    for pos in Positions:
                        if seq[pos-1] != previous_seq[pos-1]:
                            Mutation_count[pos-1] += 1
                            Variability[pos-1] += 1/yearly_seq_count
                        if seq[pos-1] == '-':
                            Gap_count[pos-1] += 1
                        elif seq[pos-1] == 'A':
                            A_count[pos-1] += 1
                        elif seq[pos-1] == 'T' or seq[pos-1] == 'U':
                            T_count[pos-1] += 1
                        elif seq[pos-1] == 'G':
                            G_count[pos-1] += 1
                        elif seq[pos-1] == 'C':
                            C_count[pos-1] += 1  
                    previous_seq = seq
    Gap_ratio = np.divide(Gap_count, total_seq_count)
    A_ratio = np.divide(A_count, total_seq_count)
    T_ratio = np.divide(T_count, total_seq_count)
    G_ratio = np.divide(G_count, total_seq_count)
    C_ratio = np.divide(C_count, total_seq_count)
    for i in range(0, seq_length):
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
                       'Mutation_count': Mutation_count,
                       'Variability': Variability,
                       'Shannon_entropy': Entropy,
                       'Gap_percentage': Gap_ratio * 100,
                       'A_percentage': A_ratio * 100,
                       'T_percentage': T_ratio * 100,
                       'G_percentage': G_ratio * 100,
                       'C_percentage': C_ratio * 100})
    df['Position'] = df['Position'].astype('int64')
    df['Mutation_count'] = df['Mutation_count'].astype('int64')
    return df

for Type in ['A','B']:
    for Segment in list('12345678'):
        globals()['complete_df_'+Type+'_'+Segment] = complete_segment_df_maker(Type, Segment)
        globals()['complete_df_'+Type+'_'+Segment].to_csv('./Data/Dataframes/Complete/complete_df_'+Type+'_'+Segment+'.csv', index = False)   


#%% Complete dataframe per type

for Type in ['A','B']:
    for Segment in list('12345678'):
        with open('./Data/Dataframes/Complete/complete_df_'+Type+'_'+Segment+'.csv') as file:
            globals()['complete_df_'+Type+'_'+Segment] = pd.read_csv(file)


def complete_type_df_maker(Type):
    df = pd.DataFrame({'Position': [],
                       'Segment': [],
                       'Segment_position': [],
                       'Mutation_count': [],
                       'Variability': [],
                       'Shannon_entropy': [],
                       'Gap_percentage': [],
                       'A_percentage': [],
                       'T_percentage': [],
                       'G_percentage': [],
                       'C_percentage': []})
    length = 0
    for Segment in list('12345678'):
        df_Segment = globals()['complete_df_'+Type+'_'+Segment].copy()
        df_Segment['Segment'] = [int(Segment)] * len(df_Segment)
        df_Segment = df_Segment.rename({'Position': 'Segment_position'}, axis='columns')
        df_Segment['Position']  = df_Segment['Segment_position'] + length
        length = max(df_Segment['Position'])
        df = df.append(df_Segment)
    df['Position'] = df['Position'].astype('int64')
    df['Segment'] = df['Segment'].astype('int64')
    df['Segment_position'] = df['Segment_position'].astype('int64')
    df['Mutation_count'] = df['Mutation_count'].astype('int64')
    return df

for Type in ['A','B']:
    globals()['complete_df_'+Type] = complete_type_df_maker(Type)
    globals()['complete_df_'+Type].to_csv('./Data/Dataframes/Complete/complete_df_'+Type+'.csv', index = False)   

#%% Ungapped dataframe per segment

def segment_df_ungapper(Type, Segment):
    df = globals()['complete_df_'+Type+'_'+Segment].copy()
    df = df[df['Gap_percentage'] < 90].sort_values('Position')
    df = df.rename({'Position': 'Aligned_position'}, axis='columns')
    df['Position'] = list(range(1, len(df)+1))
    return df

for Type in ['A', 'B']:
    for Segment in list('12345678'):
        globals()['df_'+Type+'_'+Segment] = segment_df_ungapper(Type, Segment)
        globals()['df_'+Type+'_'+Segment].to_csv('./Data/Dataframes/Ungapped/df_'+Type+'_'+Segment+'.csv', index = False)
        
#%% Ungapped dataframe per type

def type_df_ungapper(Type):
    df = globals()['complete_df_'+Type].copy()
    df = df[df['Gap_percentage'] < 90].sort_values('Position')
    df = df.rename({'Position': 'Aligned_position', 'Segment_position': 'Aligned_segment_position'}, axis='columns')
    df['Position'] = list(range(1, len(df)+1))
    Segment_position = []
    for Segment in list('12345678'):
        df_Segment = df[df['Segment'] == int(Segment)]
        Segment_position += list(range(1, len(df_Segment)+1))
    df['Segment_position'] = Segment_position
    return df

for Type in ['A', 'B']:
    globals()['df_'+Type] = type_df_ungapper(Type)
    globals()['df_'+Type].to_csv('./Data/Dataframes/Ungapped/df_'+Type+'.csv', index = False)
        