#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct  2 12:08:06 2021

@author: mennovandamme
"""

#### Clean and store the sequence data and metadata


#%% Dependencies

import re, json, os, subprocess
# conda install -c anaconda biopython
from Bio import SeqIO

#%% Directories

if not os.path.isdir('./Data'): os.mkdir('./Data')
if not os.path.isdir('./Data/NCBI'): os.mkdir('./Data/NCBI')
if not os.path.isdir('./Data/Sorted'): os.mkdir('./Data/Sorted')
if not os.path.isdir('./Data/Clean'): os.mkdir('./Data/Clean')
if not os.path.isdir('./Data/Metadata'): os.mkdir('./Data/Metadata')

#%% Downloading Influenza NCBI database using rsync

# donwload date: 11th of May

## Downloading the fasta and metadata files via Rsync
# The (GNU zipped) files we need:
#   influenza_na.dat.gz - Table with supplementary nucleotide data
#   influenza.fna.gz - FASTA nucleotide

cmd = 'rsync -aP --include=influenza_na.dat.gz --include=influenza.fna.gz --exclude=* rsync://ftp.ncbi.nlm.nih.gov/genomes/INFLUENZA/ ./Data/NCBI'
subprocess.run(cmd, shell = True)

## Unzipping the files

cmd = 'gunzip --force -r ./Data/NCBI'
subprocess.run(cmd, shell = True)

#%% Reading in the metadata

metadata = {}

## regex to extract the virus type (A or B) from the virus name
# Virus_name (column 7) looks like: 
#   A: Influenza A virus (A/Arequipa/FLU3845/2006(H3)) 
#   B: Influenza B virus (B/Jiangsu/10e9/2003)
#   problem some viruses are unidentified: unidentified influenza virus
Type_re = re.compile('^Influenza ([AB])')

## regex to check for a correct date
# some collection dates are year/month/date: 2011/01/20
# others are year/month: 2011/05
# others are : 1985
# unwanted collection dates: empty ones ('') or 'unknown'
Date_re = re.compile('^([0-9]+/*[0-9]*/*[0-9]*)')
Month_re = re.compile('^[0-9]+/([0-9]+)/*[0-9]*')
Day_re = re.compile('^[0-9]+/[0-9]+/([0-9]+)')

# list of IDs without a known type (A or B)
no_type_IDs = []

# list of IDs without a complete date
incomplete_date_IDs = []

# list of IDs that are not from human host
other_host_IDs = []

# list of IDs that aren't complete sequences
nc_IDs = []

t0 = time.time()
with open('./Data/NCBI/influenza_na.dat', 'r') as file:
    for line in file:
        record = line.replace('\n','').split('\t')  
        ID = record[0]
        if Type_re.search(record[7]):
            Type = Type_re.search(record[7]).group(1)
        else:
            Type = 'unknown'
            no_type_IDs.append(ID)
        if Date_re.search(record[5]):
            Date = Date_re.search(record[5]).group(1)
            Year = int(Date[0:4])
            if Month_re.search(record[5]):
                Month = int(Month_re.search(record[5]).group(1))
                if Day_re.search(record[5]):
                    Day = int(Day_re.search(record[5]).group(1))
                else:
                    Day = 'unknown'
                    incomplete_date_IDs.append(ID)
            else:
                 Month = 'unknown'
                 Day = 'unknown'
                 incomplete_date_IDs.append(ID)
        else:
            Date = 'unknown'
            Year = 'unknown'
            Month = 'unknown'
            Day = 'unknown'
            incomplete_date_IDs.append(ID)
        Host = record[1]
        if Host != 'Human':
            other_host_IDs.append(ID)
        Complete = record[10]
        if Complete != 'c':
            nc_IDs.append(ID)
        metadata[ID] = {'Host': Host,
                        'Type': Type,
                        'Segment': record[2],
                        'Subtype': record[3],
                        'Country': record[4],
                        'Collection_date': Date,
                        'Year': Year,
                        'Month': Month,
                        'Day': Day,
                        'Seq_length': record[6],
                        'Virus_name': record[7], 
                        'Completeness': Complete,
                        'Description': ''
                        }

len(metadata)
# 817587

len(no_type_IDs)
# 2142

len(incomplete_date_IDs)
# 132575

len(other_host_IDs)
# 324341

len(nc_IDs)
# 119520

t1 = time.time()
print('Reading in the metadata took:', round(t1-t0), 'seconds')

#%% Reading in the sequence data and sorting them into separate files per type and segment

## regex to extract the genbank acession ID from the seq_record.id/FASTA header
# FASTA header looks like this:
#   >gi|58576|gb|X52226|Influenza A virus (A/FPV/Rostock/34(H7N1)) gene for neuraminidase, genomic RNA
# sequence record ID looks like this (first part of FASTA header):
#   gi|222425131|gb|AB451147|Influenza
ID_re = re.compile('^gi\|[0-9]+\|gb\|(.*)\|.*$')

# set of unwanted IDs
bad_IDs = set(no_type_IDs + incomplete_date_IDs + other_host_IDs + nc_IDs)

def sequence_sorter(fasta_input, metadata, bad_IDs):   
    for Type in ['A', 'B']:
        for Segment in list('12345678'):
            locals()[Type+'_'+Segment] = []   
    # sorting the sequence records into the correct lists
    for seq_record in SeqIO.parse(fasta_input, 'fasta'):
        ID = ID_re.search(seq_record.id).group(1)
        seq_record.id = ID
        metadata[ID]['Description'] = seq_record.description
        seq_record.description = ''
        if ID not in bad_IDs:
            Type = metadata[seq_record.id]['Type']
            Segment = metadata[seq_record.id]['Segment']
            locals()[Type+'_'+Segment].append(seq_record)            
    # writing the sequence records into their fasta files & printing the counts
    print('Sorted sequence counts:')
    for Type in ['A', 'B']:
        for Segment in list('12345678'):
            SeqIO.write(locals()[Type+'_'+Segment], './Data/Sorted/'+Type+'_'+Segment+'.fasta', 'fasta') 
            print(Type+'_'+Segment+':', len(locals()[Type+'_'+Segment]))

sequence_sorter('./Data/NCBI/influenza.fna', metadata, bad_IDs)     

# Sorted sequence counts BEFORE REMOVING INCOMPLETE DATE IDS:
# A_1: 37111
# A_2: 36421
# A_3: 37378
# A_4: 53777
# A_5: 38357
# A_6: 47288
# A_7: 42567
# A_8: 38602
# B_1: 11104
# B_2: 11145
# B_3: 11118
# B_4: 16146
# B_5: 11220
# B_6: 12743
# B_7: 11245
# B_8: 12411

# Sorted sequence counts AFTER REMOVING INCOMPLETE DATE IDS:
# A_1: 34556
# A_2: 33826
# A_3: 34792
# A_4: 49567
# A_5: 35675
# A_6: 43561
# A_7: 39422
# A_8: 35806
# B_1: 10933
# B_2: 10973
# B_3: 10943
# B_4: 15398
# B_5: 11037
# B_6: 12378
# B_7: 11043
# B_8: 12188              

with open('./Data/Metadata/metadata_complete.json', 'w') as file:
    json.dump(metadata, file)

#%% Removing ambiguous IDs

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

# list of IDs of sequences with ambiguous characters
ambi_IDs = []

# dict that keeps track of ambiguous sequence counts
ambi_counts = {}
for Type in ['A', 'B']:
    ambi_counts[Type] = {}
    for Segment in list('12345678'):
        ambi_counts[Type][Segment] = 0
        
# list of IDs of duplicate sequences
duplicate_IDs = []

# dict that keeps track of which sequence (oldest ID) has which duplicates (other IDs)
dupes = {}

# dict that keeps track of duplicate sequence counts
dupe_counts = {}
for Type in ['A', 'B']:
    dupe_counts[Type] = {}
    for Segment in list('12345678'):
        dupe_counts[Type][Segment] = 0

def sequence_cleaner(Type, Segment, metadata):
    '''
    removes sequences which have too many ambiguous characters
    removes duplicate sequences
    writes a new fasta file ontaining the cleaned up sequences with the oldest ID
    '''
    sequences = {}
    input_file = './Data/Sorted/'+Type+'_'+Segment+'.fasta'
    output_file = './Data/Clean/'+Type+'_'+Segment+'.fasta'    
    for seq_record in SeqIO.parse(input_file, 'fasta'):
        ID = seq_record.id
        seq = str(seq_record.seq).upper()        
        # if any charcter in the sequence is ambiguous append to list & count the amount per type/segment
        if any([i not in ['A', 'T', 'G', 'C'] for i in set(seq)]):
            ambi_IDs.append(ID)
            ambi_counts[Type][Segment] += 1
        else:
            # add sequence as dict key (if not in dict) and the ID to value list
            if seq not in sequences: 
                sequences[seq] = [ID]
            else:
                sequences[seq].append(ID)
                dupe_counts[Type][Segment] += 1   
    print(Type+'_'+Segment+':', len(sequences))                                   
    with open(output_file, 'w') as f:
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
            f.write('>' + oldest_ID + '\n' + seq + '\n')    
    

print('Cleaned sequence counts:')
for Type in ['A', 'B']:
    for Segment in list('12345678'):
        sequence_cleaner(Type, Segment, metadata)

# Cleaned sequence counts:
# A_1: 20188
# A_2: 19603
# A_3: 19793
# A_4: 29435
# A_5: 16980
# A_6: 23232
# A_7: 13121
# A_8: 13133
# B_1: 6320
# B_2: 6162
# B_3: 6446
# B_4: 8776
# B_5: 5874
# B_6: 6778
# B_7: 4540
# B_8: 5559

# ambi_counts
# {'A': 
#  {'1': 2011,
#   '2': 2016,
#   '3': 1999,
#   '4': 3440,
#   '5': 1269,
#   '6': 2892,
#   '7': 1084,
#   '8': 1085},
#  'B': 
#  {'1': 594,
#   '2': 572,
#   '3': 582,
#   '4': 684,
#   '5': 520,
#   '6': 544,
#   '7': 280,
#   '8': 338}}

# dupe_counts
# {'A': 
#  {'1': 12357,
#   '2': 12207,
#   '3': 13000,
#   '4': 16692,
#   '5': 17426,
#   '6': 17437,
#   '7': 25217,
#   '8': 21588},
#  'B': 
#  {'1': 4019,
#   '2': 4239,
#   '3': 3915,
#   '4': 5938,
#   '5': 4643,
#   '6': 5056,
#   '7': 6223,
#   '8': 6291}}

#%% Making and storing a clean metadata dictionary

bad_IDs = set(no_type_IDs + incomplete_date_IDs + other_host_IDs + nc_IDs + duplicate_IDs + ambi_IDs)

metadata_clean = metadata.copy()

for ID in bad_IDs:
    del metadata_clean[ID]

for ID in metadata_clean:
    metadata_clean[ID]['Duplicate_IDs'] = dupes[ID]
    metadata_clean[ID]['Duplicate_count'] = len(dupes[ID])
  
len(metadata_clean)
# 376098
    
with open('./Data/Metadata/metadata_clean.json', 'w') as file:
    json.dump(metadata_clean, file, default = str)

# restructuring metadata
metadata_new = {}
for Type in ['A', 'B']:
    metadata_new[Type] = {}
    for Segment in list('12345678'):
        metadata_new[Type][Segment] = {}

for ID in metadata_clean:
    Type = metadata_clean[ID]['Type']
    Segment = metadata_clean[ID]['Segment']
    metadata_new[Type][Segment][ID] = metadata_clean[ID]
    
with open('./Data/Metadata/metadata.json', 'w') as file:
    json.dump(metadata_new, file, default = str)