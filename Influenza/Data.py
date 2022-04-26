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
if not os.path.isdir('./Data/NCBI'): os.mkdir('./Data/NCBI')
if not os.path.isdir('./Data/Sorted'): os.mkdir('./Data/Sorted')
if not os.path.isdir('./Data/Clean'): os.mkdir('./Data/Clean')
if not os.path.isdir('./Data/Metadata'): os.mkdir('./Data/Metadata')

#%% Regexes

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

## regex to extract the genbank acession ID from the seq_record.id/FASTA header
# FASTA header looks like this:
#   >gi|58576|gb|X52226|Influenza A virus (A/FPV/Rostock/34(H7N1)) gene for neuraminidase, genomic RNA
# sequence record ID looks like this (first part of FASTA header):
#   gi|222425131|gb|AB451147|Influenza

ID_re = re.compile('^gi\|[0-9]+\|gb\|(.*)\|.*$')

#%% Downloading Influenza NCBI database using rsync

t0 = time.time()

## Downloading the fasta and metadata files via Rsync
# The (GNU zipped) files we need:
#   influenza_na.dat.gz - Table with supplementary nucleotide data
#   influenza.fna.gz - FASTA nucleotide

cmd = 'rsync -aP --include=influenza_na.dat.gz --include=influenza.fna.gz --exclude=* rsync://ftp.ncbi.nlm.nih.gov/genomes/INFLUENZA/ ./Data/NCBI'
subprocess.run(cmd, shell = True)

## Unzipping the files

cmd = 'gunzip --force -r ./Data/NCBI'
subprocess.run(cmd, shell = True)

t1 = time.time()
print('Downloading Influenza NCBI database took', round((t1-t0)/60, 1), 'minutes')
# 

#%% Reading in the metadata

## loading in the metadata from influenza_na.dat
# creating a Type and Year key
# if the Type or Date is doesn't get recognize by regex set to 'unknown'

metadata = {}
no_date_IDs = []
no_type_IDs = []
other_host_IDs = []
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
            else:
                 Month = 'unknown'
                 Day = 'unknown'
        else:
            Date = 'unknown'
            Year = 'unknown'
            Month = 'unknown'
            Day = 'unknown'
            no_date_IDs.append(ID)
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

t1 = time.time()
print('Reading in the metadata took:', round(t1-t0), 'seconds')

#%% Reading in the sequence data and sorting them into separate files

## make a list with the IDs that don't have an identified type or collection date
bad_IDs = set(no_date_IDs + no_type_IDs + other_host_IDs + nc_IDs)

## sequence_sorter function:
#   gives the sequences a new ID that is only the GenBank accession ID
#   stores the sequence record description in the metadata dictionary
#   remove the sequence record description so it doesn't get put in the output FASTA header
#   sort them by influenza type (A or B) and by segment (1-8)
#   write to separate fasta files to an output directory
#   print the number of sequences in each file

def sequence_sorter(fasta_input, metadata):
    for Type in ['A', 'B']:
        for Segment in list('12345678'):
            locals()[Type+'_'+Segment] = []
    
    # sorting the sequence records into the correct lists
    for seq_record in SeqIO.parse(fasta_input, 'fasta'):
        ID = ID_re.search(seq_record.id).group(1)
        seq_record.id = ID
        if ID not in bad_IDs:
            metadata[seq_record.id]['Description'] = seq_record.description
            seq_record.description = ''
            Host = metadata[seq_record.id]['Host']
            Type = metadata[seq_record.id]['Type']
            Segment = metadata[seq_record.id]['Segment']
            locals()[Type+'_'+Segment].append(seq_record)
            
    # writing the sequence records into their fasta files   
    for Type in ['A', 'B']:
        for Segment in list('12345678'):
            SeqIO.write(locals()[Type+'_'+Segment], './Data/Sorted/'+Type+'_'+Segment+'.fasta', 'fasta')
    
    # print the number of sequences in each file
    for Type in ['A', 'B']:
        for Segment in list('12345678'):
            print('\t'+Host+'_'+Type+'_'+Segment+':', len(locals()[Type+'_'+Segment]))

print('Sorted sequence counts:')
sequence_sorter(fasta_input = './Data/NCBI/influenza.fna', metadata = metadata)                    

with open('./Data/Metadata/metadata_complete.json', 'w') as file:
    json.dump(metadata, file)

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

ambi_IDs = []
ambi_counts = {}
for Type in ['A', 'B']:
    ambi_counts[Type] = {}
    for Segment in list('12345678'):
        ambi_counts[Type][Segment] = 0
duplicate_IDs = []
dupes = {}

## sequence_cleaner function:
#   removes sequences which have too many ambiguous characters
#   removes duplicate sequences
#   returns a dictionary whith the oldest ID as key and it's duplicate IDs as values
#   writes a new fasta file ontaining the cleaned up sequences with the oldest ID
#   prints amount of sequences in the cleaned file

def sequence_cleaner(Type, Segment, metadata):
    sequences = {}
    input_file = './Data/Sorted/'+Type+'_'+Segment+'.fasta'
    output_file = './Data/Clean/'+Type+'_'+Segment+'.fasta'
    
    for seq_record in SeqIO.parse(input_file, 'fasta'):
        ID = seq_record.ID
        seq = str(seq_record.seq).upper()
        
        # if any charcter in the sequence is ambiguous append to list & count the amount per type/segment
        if any([i not in ['A', 'T', 'G', 'C'] for i in set(seq)]):
            ambi_IDs.append(ID)
            ambi_counts[Type][Segment] += 1
        else:
            # add sequence as dict key (if not in dict) and the ID to value list
            if seq not in sequences: 
                sequences[seq] = []
            sequences[seq].append(ID) 
            
    with open(output_file, 'w') as output_file:
        for seq in sequences:
            IDs = sequences[seq]
            if len(IDs) == 1:
                oldest_ID = IDs[0]
                dupes[oldest_ID] = 
            else:
                temp = {}
                for ID in IDs:
                    temp[ID] = metadata[ID]['Year']
                oldest_ID = min(temp.items(), key = lambda item: item[1])[0]
                duplicates = IDs.copy()
                duplicates.remove(oldest_ID)
                duplicate_IDs.extend(duplicates)
                dupes[oldest_ID] = duplicates            
            output_file.write('>' + oldest_ID + '\n' + seq + '\n')        
    print('\t'+Host+'_'+Type+'_'+Segment+':', len(sequences))

t0 = time.time()
print('Cleaned sequence counts:')
for Host in ['Human', 'Other']:
    for Type in ['A', 'B']:
        for Segment in list('12345678'):
                sequence_cleaner(Host = Host, Type = Type, Segment = Segment,
                                 metadata = metadata,
                                 min_length = min_lengths[Type][Segment], 
                                 max_ambiguous = 0)
t1 = time.time()
print('Cleaning up the sequences took:', round(t1-t0), 'seconds')


#%% Making and storing a clean metadata dictionary

t0 = time.time()

## remove the IDs which:
#   didn't meet the sequence_sorter requirements:have unidentified type or no collection date
#   didn't meet the sequence_cleaner requirements: have Ns or are shorter than min_length

bad_IDs = bad_IDs_ss + bad_IDs_sc + duplicate_IDs
metadata_clean = metadata.copy()
for ID in bad_IDs:
    del metadata_clean[ID]

for ID in metadata_clean:
    metadata_clean[ID]['Duplicate_IDs'] = dupes[ID]
    metadata_clean[ID]['Duplicate_count'] = len(dupes[ID])
    

## writing the dictionary to a .JSON file
    
with open('./Data/Metadata/metadata_clean.json', 'w') as file:
    json.dump(metadata_clean, file, default = str)

t1 = time.time()
print('Making and storing a clean metadata dictionary took:', round(t1-t0), 'seconds')

t3 = time.time()
print('Data.py took', round((t3-t2)/60, 1), 'minutes')
