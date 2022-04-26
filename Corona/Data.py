#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct  2 12:08:06 2021

@author: mennovandamme
"""

#### Clean and store the sequence data and metadata


#%% Dependencies

import re, json, time, os, subprocess, math, pickle
# conda install -c anaconda biopython
from Bio import SeqIO, Entrez

#%% Directories

if not os.path.isdir('./Data'): os.mkdir('./Data')
if not os.path.isdir('./Data/Sorted'): os.mkdir('./Data/Sorted')
if not os.path.isdir('./Data/Clean'): os.mkdir('./Data/Clean')
if not os.path.isdir('./Data/Metadata'): os.mkdir('./Data/Metadata')

#%% Regexes



#%% Downloading SARS-CoV-2 NCBI database using datasets

cmd = 'datasets download virus genome taxon SARS-CoV-2 --host human --complete-only --exclude-cds --exclude-protein --filename ./Data/NCBI.zip'
subprocess.run(cmd, shell = True)

cmd = 'unzip -o ./Data/NCBI.zip'
subprocess.run(cmd, shell = True)

#%% Reading in the metadata

## regex to check for a correct date
# some collection dates are year-month-date: 2011-01-20
# others are year-month: 2011-05
# others are year: 1985
# unwanted collection dates: empty ones ('') or 'unknown' or 'missing'

Date_re = re.compile('^([0-9]+-*[0-9]*-*[0-9]*)')
Month_re = re.compile('^[0-9]+-([0-9]+)-*[0-9]*')
Day_re = re.compile('^[0-9]+-[0-9]+-([0-9]+)')

metadata = {}

no_date_IDs = []


with open('./Data/NCBI/', 'r') as file:
    for line in file:
        record = line.replace('\n','').split('\t')  
        ID = record[0]
        if Date_re.search(record[5]):
            Date = Date_re.search(record[5]).group(1)
            Year = int(Date[0:4])
            if Month_re.search(record[5]):
                Month = int(Month_re.search(record[5]).group(1))
                if Day_re.search(record[5]):
                    Day = int(Day_re.search(record[5]).group(1))
                else:
                    Day = 'unknown'
            else:
                 Month = 'unknown'
                 Day = 'unknown'
        else:
            Date = 'unknown'
            Year = 'unknown'
            Month = 'unknown'
            Day = 'unknown'
            no_date_IDs.append(ID)
        metadata[ID] = {'Subtype': record[3],
                        'Country': record[4],
                        'Collection_date': Date,
                        'Year': Year,
                        'Month': Month,
                        'Day': Day,
                        'Seq_length': record[6],
                        'Virus_name': record[7], 
                        'Description': ''
                        }

print('Amount of sequences without date:', len(no_date_IDs))
# 
#%% Cleaning up the fasta files: 0 Ns and min length

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
            
            # if the day is unknown set to middle of the month (15th)
            # if month is unknown set to middle of the year (1st of July)
            if Month == 'unknown':
                if Day == 'unknown':
                    Month = 7
                    Day = 1
                else:
                    Day = 15

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
                    IDs_sorted.append(ID)  
    return IDs_sorted


## sequence_cleaner function:
#   removes sequences with ambiguous characters
#   removes duplicate sequences
#   returns a dictionary whith the oldest ID as key and it's duplicate IDs as values
#   writes a new fasta file ontaining the cleaned up sequences with the oldest ID
#   prints amount of sequences in the cleaned file

def sequence_cleaner(in_file, out_file, metadata):
    ambi_IDs = []
    duplicate_IDs = []
    dupes = {}
    # dictionary to hold unique sequences as keys and IDs as values
    sequences = {}
    for seq_record in SeqIO.parse(in_file, 'fasta'):
        ID = seq_record.ID
        seq = str(seq_record.seq).upper()
        metadata[ID]['Description'] = seq_record.description
        # if any charcter in the sequence is ambiguous append to list & count the amount per type/segment
        if any([i not in ['A', 'T', 'G', 'C'] for i in set(seq)]):
            ambi_IDs.append(ID)
        else:
            # add sequence as dict key (if not in dict) and the ID to value list
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
                # find the oldest ID
                IDs_sorted = sort_IDs_by_date(IDs, metadata)
                oldest_ID = IDs_sorted[0]
                duplicates = IDs_sorted[1:-1]
                duplicate_IDs.extend(duplicates)
                dupes[oldest_ID] = duplicates            
            output_file.write('>' + oldest_ID + '\n' + seq + '\n')  
    return ambi_IDs, duplicate_IDs, dupes
            

ambi_IDs, duplicate_IDs, dupes = sequence_cleaner('./Data/NCBI/', './Data/NCBI/Clean/SARSCoV2.fasta', metadata)

len(ambi_IDs)
# 

len(duplicate_IDs)
# 

#%% Making a clean metadata dictionary

bad_IDs = no_date_IDs + ambi_IDs + duplicate_IDs
metadata_clean = metadata.copy()

for ID in bad_IDs:
    del metadata_clean[ID]
    
len(metadata_clean)
# 

for ID in metadata_clean:
    metadata_clean[ID]['Duplicate_IDs'] = dupes[ID]
    metadata_clean[ID]['Duplicate_count'] = len(dupes[ID])
    
with open('./Data/Metadata/complete_metadata.json', 'w') as file:
    json.dump(metadata, file, default = str)

with open('./Data/Metadata/metadata.json', 'w') as file:
    json.dump(metadata_clean, file, default = str)
