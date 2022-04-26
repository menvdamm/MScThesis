#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct  2 12:14:06 2021

@author: mennovandamme
"""

#### Restructure clean metadata and find the reference sequences for the mutation plot

#%% Dependencies

from Bio import SeqIO
import json, re

#%% Loading the necessary data
    
with open('./Data/Metadata/metadata_clean.json') as file:
    metadata_clean = json.load(file)

with open('./Data/Metadata/references.json') as file:
    references = json.load(file)
    
#%% Restructuring metadata for human host (and storing in the aligned reference sequencees)
metadata_human = {}

for Type in ['A', 'B']:
    metadata_human[Type] = {}
    for Segment in list('12345678'):
        metadata_human[Type][Segment] = {}
        file = './Data/Aligned/Human/'+Type+'_'+Segment+'.fasta'
        for seq_record in SeqIO.parse(file, 'fasta'):
            metadata_human[Type][Segment][seq_record.id] = metadata_clean[seq_record.id]
            metadata_human[Type][Segment][seq_record.id]['Aligned_sequence'] = seq_record.seq
            if seq_record.id == references[Type][Segment]['Alt_ID']:
                references[Type][Segment]['Aligned_sequence'] = seq_record.seq

#%% Storing the data

## writing the dictionaries into .JSON files
    
with open('./Data/Metadata/metadata_human.json', 'w') as file:
    json.dump(metadata_human, file, default = str)

with open('./Data/Metadata/references.json', 'w') as file:
    json.dump(references, file, default=str)