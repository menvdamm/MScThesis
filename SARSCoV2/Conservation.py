#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 15 17:54:48 2021

@author: mennovandamme
"""

### Creating a dataframe for all data

#%% Dependencies

import json, os, subprocess
import pandas as pd
import numpy as np
from Bio import SeqIO
from scipy import stats

#%% Directories

if not os.path.isdir('./Data/Dataframes'): os.mkdir('./Data/Dataframes')
if not os.path.isdir('./Data/Genome'): os.mkdir('./Data/Genome')

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
        
#%% Complete dataframe

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

#%% Ungapped dataframe

def ungap_df(df, allowed_gap_percentage):
    df = df[df['Gap_percentage'] < allowed_gap_percentage].sort_values('Position')
    df['Position'] = list(range(1, len(df)+1))
    return df

df = ungap_df(complete_df, allowed_gap_percentage)
df.to_csv('./Data/Dataframes/df.csv', index = False)

#%% Scoring 30 nucleotide regions

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
    score_df = score_df.sort_values('Begin_position')
    return score_df  

complete_score_df = score_window(df, 30)
complete_score_df.to_csv('./Data/Dataframes/complete_score_df.csv', index = False)

#%% Making consensus sequence 

with open('./Data/Genome/SARSCoV2_consensus.fasta', 'w') as file:
    consensus = ''
    for row in df.sort_values('Position').iterrows():
        percentages = {'A': row[1]['A_percentage'],
                       'T': row[1]['T_percentage'],
                       'G': row[1]['G_percentage'],
                       'C': row[1]['C_percentage']}
        most_occuring = max(percentages, key = percentages.get)
        consensus += most_occuring
    file.write('>SARS-CoV-2' + '\n' + consensus)

#%% Align CDS to consensus

# manually downloaded coding sequences from NCBI SARS-CoV-2 reference sequence
# https://www.ncbi.nlm.nih.gov/nuccore/1798174254?report=genbank
# => './Data/Genome/CDS.fasta'

CDS_file = './Data/Genome/CDS.fasta'
consensus_file = './Data/Genome/SARSCoV2_consensus.fasta'
output_file = './Data/Genome/SARSCoV2.fasta'
cmd = 'mafft --auto  --preservecase --keeplength --addfragments '+CDS_file+' --reorder --thread -1 '+consensus_file+' > '+output_file
subprocess.run(cmd, shell = True)

#%% Extract positions of coding sequences

begin_pos = []
end_pos = []
for seq_record in SeqIO.parse('./Data/Genome/SARSCoV2.fasta', 'fasta'):
    ID = seq_record.id
    if ID != 'SARS-CoV-2':
        seq = str(seq_record.seq)
        for i in range(len(seq)):
             if seq[i] != '-':
                 begin_pos.append(i+1)
                 break
        for i in range(len(seq)-1,0,-1):
             if seq[i] != '-':
                 end_pos.append(i+1)
                 break 
                
CDS_df = pd.DataFrame({'Begin_position': begin_pos,
                       'End_position': end_pos,
                       'Protein': ['ORF1ab', 'ORF1a', 'S', '3a', 'E', 'M', '6', '7a', '7b', '8', 'N', '10']})            

CDS_df.to_csv('./Data/Genome/CDS_df.csv', index = False)

#%% Filter out any 30bp regions that aren't completely in the coding regions

def filter_scores(df, CDS_df):
    indices = []
    proteins = []
    for row in df.iterrows():
        b = int(row[1]['Begin_position'])
        e = int(row[1]['End_position'])  
        for i in range(len(CDS_df)):
            CDSb = CDS_df['Begin_position'][i]
            CDSe = CDS_df['End_position'][i]
            if CDSb <= b and e <= CDSe:
                proteins.append(CDS_df['Protein'][i])
                indices.append(row[0])
                break
    score_df = df.loc[indices]
    score_df['Protein'] = proteins
    return score_df  

score_df = filter_scores(complete_score_df, CDS_df)
score_df.to_csv('./Data/Dataframes/score_df.csv', index = False)

#%% make csv for the top ten sequences per metric

for metric in ['Shannon_entropy', 'Mutability']:
    sequences = []
    small_score_df = score_df.sort_values(metric)[0:10]
    with open('./Data/Genome/SARSCoV2_consensus.fasta', 'r') as f:
        consensus = f.readlines()[1]
        for row in small_score_df.iterrows():
            b = int(row[1]['Begin_position'])
            e = int(row[1]['End_position'])  
            seq = consensus[b:e+1]
            sequences.append(seq)
    small_score_df['Positions'] = [str(b) + '-' + str(e) for b,e in zip(small_score_df['Begin_position'], small_score_df['End_position'])]
    small_score_df['Sequence'] = sequences
    small_score_df['Shannon_entropy'] = round(small_score_df['Shannon_entropy'], 5)
    small_score_df['Mutability'] = round(small_score_df['Mutability'], 5)
    small_score_df['Gap_percentage'] = round(small_score_df['Gap_percentage'], 5)
    small_score_df.to_csv('./Data/Dataframes/small_score_df_'+metric+'.csv', index = False) 

            
#%% average shannon entropy            
            
print(sum(df['Shannon_entropy'])/len(df['Shannon_entropy']))
 
print(sum(df['Mutability'])/len(df['Mutability']))      

#%% average & median Shannon entropy & mutability per CDS

# adding protein to conservation df per position
proteins = []
for row in df.iterrows():
    pos = row[1]['Position']
    CDS = False
    for i in range(len(CDS_df)):
        CDSb = CDS_df['Begin_position'][i]
        CDSe = CDS_df['End_position'][i]
        if CDSb <= pos and pos <= CDSe:
            proteins.append(CDS_df['Protein'][i])
            CDS = True
            break
    if CDS == False:
        proteins.append('UTR')
df['Protein'] = proteins
df.to_csv('./Data/Dataframes/df.csv', index = False)

for protein in list(set(df['Protein'])):
    print(protein)
    print('\tAverage Shannon entropy: ', round(np.average(df['Shannon_entropy'][df['Protein'] == protein]), 5))
    print('\tAverage Mutability: ', round(np.average(df['Mutability'][df['Protein'] == protein]), 5))
    print('\tMedian Shannon entropy: ', round(np.median(df['Shannon_entropy'][df['Protein'] == protein]), 5))
    print('\tMedian Mutability: ', round(np.median(df['Mutability'][df['Protein'] == protein]), 5))

print('Entire genome')
print('\tAverage Shannon entropy: ', round(np.average(df['Shannon_entropy']), 5))
print('\tAverage Mutability: ', round(np.average(df['Mutability']), 5))
print('\tMedian Shannon entropy: ', round(np.median(df['Shannon_entropy']), 5))
print('\tMedian Mutability: ', round(np.median(df['Mutability']), 5)) 
   
#%% correlation between Shannon entropy and mutability

print(stats.spearmanr(df['Shannon_entropy'], df['Mutability']))

print(stats.spearmanr(score_df['Shannon_entropy'], score_df['Mutability']))


