#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct  2 12:08:06 2021

@author: mennovandamme
"""

#### Download, clean and store the sequence data and metadata


#%% Dependencies

import re, json, os, subprocess
# conda install -c anaconda biopython
from Bio import SeqIO

#%% Directories

if not os.path.isdir('./Data'): os.mkdir('./Data')
if not os.path.isdir('./Data/NCBI'): os.mkdir('./Data/NCBI')
if not os.path.isdir('./Data/Clean'): os.mkdir('./Data/Clean')
if not os.path.isdir('./Data/Metadata'): os.mkdir('./Data/Metadata')

#%% Downloading SARS-CoV-2 NCBI database using datasets

cmd = 'datasets download virus genome taxon SARS-CoV-2 --host human --complete-only --exclude-cds --exclude-protein --filename ./Data/NCBI/ncbi_dataset.zip'
subprocess.run(cmd, shell = True)

cmd = 'unzip -o ./Data/NCBI/ncbi_dataset.zip -d ./Data/NCBI'
subprocess.run(cmd, shell = True)

#%% Reading in the metadata

# regex to check for full date
Date_re = re.compile('^[0-9]+-[0-9]+-[0-9]+')

# dictionary for the metada
metadata = {}

# list of IDs without full date
no_date_IDs = []

with open('./Data/NCBI/ncbi_dataset/data/data_report.jsonl', 'r') as file:
    for line in file:
        record = json.loads(line)
        ID = record['accession']        
        if 'isolate' in record.keys() and 'collectionDate' in record['isolate'].keys() and Date_re.search(record['isolate']['collectionDate']):
            Date = record['isolate']['collectionDate']
            ymd = Date.split('-')
            Year = int(ymd[0])
            Month = int(ymd[1])
            Day = int(ymd[2])
        else:
            Date = 'unknown'
            Year = 'unknown'
            Month = 'unknown'
            Day = 'unknown'
            no_date_IDs.append(ID)           
        if 'location' in record.keys():
            Country = record['location']['geographicLocation']
        else:
            Country = 'unknown'
        metadata[ID] = {'Country': Country,
                        'Collection_date': Date,
                        'Year': Year,
                        'Month': Month,
                        'Day': Day,
                        'Seq_length': record['length'],
                        'Virus_name': record['virus']['sciName'], 
                        'Pangolin_class': record['virus']['pangolinClassification'],
                        'Description': ''}

len(no_date_IDs)
# 27420

len(metadata)
# 1074845

#%% Removing ambiguous sequences and sequences without collection data

def remove_ambiguous_and_no_date(in_file, out_file, metadata, no_date_IDs):
    ambi_IDs = []
    with open(out_file, 'w') as output_file:
        for seq_record in SeqIO.parse(in_file, 'fasta'):
            ID = seq_record.id
            metadata[ID]['Description'] = seq_record.description
            if ID not in no_date_IDs:
                seq = str(seq_record.seq).upper()
                # if any charcter in the sequence is ambiguous append to list
                if any([i not in ['A', 'T', 'G', 'C'] for i in set(seq)]):
                    ambi_IDs.append(ID)
                else:
                    output_file.write('>' + ID + '\n' + seq + '\n')  
    return ambi_IDs

ambi_IDs = remove_ambiguous_and_no_date('./Data/NCBI/ncbi_dataset/data/genomic.fna', './Data/NCBI/SARSCoV2.fasta', metadata, no_date_IDs)

with open('./Data/Metadata/complete_metadata.json', 'w') as file:
    json.dump(metadata, file)

len(ambi_IDs)
# 580976

#%% Removing duplicate sequences

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

def remove_duplicates(in_file, out_file, metadata):
    duplicate_IDs = []
    dupes = {}
    # dictionary to hold unique sequences as keys and IDs as values
    sequences = {}
    for seq_record in SeqIO.parse(in_file, 'fasta'):
        ID = seq_record.id
        seq = str(seq_record.seq).upper()
        if seq not in sequences: 
            sequences[seq] = []
        sequences[seq].append(ID)     
    with open(out_file, 'w') as output_file:
        for seq in sequences:
            IDs = sequences[seq]
            if len(IDs) == 1:
                oldest_ID = IDs[0]
                dupes[oldest_ID] = []
            else:
                # find the oldest ID and keep that one
                IDs_sorted = sort_IDs_by_date(IDs, metadata)
                oldest_ID = IDs_sorted[0]
                duplicates = IDs_sorted[1:]
                duplicate_IDs.extend(duplicates)
                dupes[oldest_ID] = duplicates            
            output_file.write('>' + oldest_ID + '\n' + seq + '\n')  
    return duplicate_IDs, dupes
            
duplicate_IDs, dupes = remove_duplicates('./Data/NCBI/SARSCoV2.fasta', './Data/Clean/SARSCoV2.fasta', metadata)

len(duplicate_IDs)
# 117771

#%% Making a clean metadata dictionary

metadata_clean = metadata.copy()

bad_IDs = no_date_IDs + ambi_IDs + duplicate_IDs

for ID in bad_IDs:
    del metadata_clean[ID]
    
for ID in metadata_clean:
    metadata_clean[ID]['Duplicate_IDs'] = dupes[ID]
    metadata_clean[ID]['Duplicate_count'] = len(dupes[ID])

with open('./Data/Metadata/metadata.json', 'w') as f:
    json.dump(metadata_clean, f)
