#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 30 10:26:45 2022

@author: mennovandamme
"""

### Secondary RNA structure prediction

#%% Dependencies

import subprocess, os
import pandas as pd
import numpy as np

#%% Directories

if not os.path.isdir('./Data/Structure'): os.mkdir('./Data/Structure')

#%% Loading the necessary data

with open('./Data/Dataframes/df.csv', 'r') as f:
    df = pd.read_csv(f)  
    
#%% Making consensus sequence 

with open('./Data/Structure/SARSCoV2_consensus.fasta', 'w+') as file:
    seq = ''
    for row in df.sort_values('Position').iterrows():
        percentages = {'A': row[1]['A_percentage'],
                       'U': row[1]['T_percentage'],
                       'G': row[1]['G_percentage'],
                       'C': row[1]['C_percentage']}
        most_occuring = max(percentages, key = percentages.get)
        seq += most_occuring
    file.write('>SARS-CoV-2' + '\n' + seq + '\n')

#%% RNA secondary structure prediction

# Done one HPC: RNAfold.sh

# Extracting the centroid structure from the .ifold file

with open('./Data/Structure/SARSCoV2.ifold', 'r') as f:
    lines = f.readlines()
    sequence = lines[1].strip()
    structure = lines[2].strip().split(' ')[0]
            
# Making a file with sequence and structure

with open('./Data/Structure/SARSCoV2_structure.txt', 'w') as file:         
    file.write('>SARS-CoV-2' + '\n' + sequence + '\n' + structure + '\n')

#%% Making a color file for visualization in ViennaRNA webbrowser

colors = df.sort_values('Position')['Mutability'].copy()
log_colors = np.log10(colors + 1)/max(np.log10(colors + 1))
with open('./Data/Structure/colors.txt', 'w') as file:
    file.write('>SARS-CoV-2' + '\n')
    for color in log_colors:
        file.write(str(round(color, 2))+'\n')
        
#%% Transforming dotbracket structure into element structure

# using rnaConvert.py script from the ViennaRNA forgi library
# https://viennarna.github.io/forgi/scripts_tutorial.html
# https://viennarna.github.io/forgi/graph_tutorial.html
#
# element structure:
#   f = fiveprime end, the unpaired nucleotides at the 5’ end of a molecule/ chain
#   t = threeprime end, the unpaired nucleotides at the 3’ end of a molecule/ chain
#   i = interior loop, bulged out nucleotides and interior loops. An interior loop can contain unpaired bases on either strand or on both strands, flanked by stems on either side.
#   s = stem, regions of contiguous canonical Watson-Crick base-paired nucleotides. By default, stems have at least 2 consecutive basepairs.
#   m = multiloop segment, single-stranded regions bewteen two stems. In the current version of forgi, pseudo-knots and exterior loops segments between stems are treated as multiloop segments.
#   h = hairpin loop, closed loop at the end of a stem

cmd = 'rnaConvert.py ./Data/Structure/SARSCoV2_structure.txt -T element_string --filename ./Data/Structure/SARSCoV2'
subprocess.run(cmd, shell = True)

#%% Adding structure element to conservation dataframe

with open('./Data/Structure/SARSCoV2001.element_string') as file:
    df['Element'] = list(file.readlines()[1].strip())
    
    
df['Element'] = df['Element'].replace({"f": "5' unpaired",                
                                       "i": "interior loop & bulge",
                                       "s": "stem",
                                       "m": "multiloop",
                                       "h": "hairpin loop",
                                       "t": "3' unpaired"
                                       })

df.to_csv('./Data/Dataframes/df.csv', index = False)

#%% RUnning ViennaRNA forna server





                