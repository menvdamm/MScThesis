#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  7 16:46:02 2021

@author: mennovandamme
"""

#### Pipeline to create 2D folding figure of the RNAs

#%% Dependencies

import json, subprocess, shutil, re, webbrowser, os
import pandas as pd
import numpy as np
#conda install -c bioconda viennarna
#import RNA

#%% Directories

if not os.path.isdir('./Data/Structures'): os.mkdir('./Data/Structures')
if not os.path.isdir('./Data/Structures/Colors'): os.mkdir('./Data/Structures/Colors')

#%% Loading the necessary data

with open('./Data/Metadata/references.json') as file:
    references = json.load(file)

for Type in ['A', 'B']:
    with open('./Data/Dataframes/Ungapped/df_'+Type+'.csv') as file:
        globals()['df_'+Type] = pd.read_csv(file)    
    for Segment in list('12345678'): 
        with open('./Data/Dataframes/Ungapped/df_'+Type+'_'+Segment+'.csv') as file:
            globals()['df_'+Type+'_'+Segment] = pd.read_csv(file)

#%% Making ungapped synthetic reference fasta

## writing a fasta file for the ungapped reference sequences
with open('./Data/Structures/references.fasta', 'w+') as file:
    for Type in ['A', 'B']:
        for Segment in list('12345678'):
            df = globals()['df_'+Type+'_'+Segment]
            ID = references[Type][Segment]['ID']
            al_seq = references[Type][Segment]['Aligned_sequence']
            seq = ''
            for row in df.sort_values('Position').iterrows():
                pos = int(row[1]['Position'])
                if al_seq[pos-1] != '-':
                    seq += al_seq[pos-1]
                else:
                    percentages = {'A': row[1]['A_percentage'],
                                   'U': row[1]['T_percentage'],
                                   'G': row[1]['G_percentage'],
                                   'C': row[1]['C_percentage']}
                    most_occuring = max(percentages, key = percentages.get)
                    seq += most_occuring
            references[Type][Segment]['Ungapped_sequence'] = seq
            file.write('>' + ID + '_' + Type + '_' + Segment + '\n' + seq + '\n')

with open('./Data/Metadata/references.json', 'w') as file:
    json.dump(references, file, default = str)
        
#%% RNA secondary structure prediction of ungapped reference sequences

# Global RNa structure prediction:
# https://www.tbi.univie.ac.at/RNA/RNAfold.1.html
# suggested by webserver:
#
#   −−noLP  Produce structures without lonely pairs (helices of length 1).
#           (default=off) For partition function folding this only disallows 
#           pairs that can only occur isolated. Other pairs may still
#           occasionally occur as helices of length 1.    
#       
#   -p      Calculate the partition function and base pairing probability matrix.
#           This notation makes use of the letters " . , | { } ( ) ".
#           On the next line the centroid structure as derived from the pair 
#           probabilities together with its free energy and distance to the ensemble is shown.
#           Note that unless you also specify −d2 or −d0, the partition function 
#           and mfe calculations will use a slightly different energy model.
#           An additionally passed value to this option changes the behavior of 
#           partition function calculation: −p0 Calculate the partition function 
#           but not the pair probabilities, saving about 50% in runtime.
#
#   -d      How to treat "dangling end" energies for bases adjacent to helices 
#           in free ends and multi−loops (default=‘2’). With −d1 only unpaired 
#           bases can participate in at most one dangling end. With −d2 this 
#           check is ignored, dangling energies will be added for the bases 
#           adjacent to a helix on both sides in any case; this is the default 
#           for mfe and partition function folding (−p). The option −d0 ignores 
#           dangling ends altogether (mostly for debugging). With −d3 mfe folding 
#           will allow coaxial stacking of adjacent helices in multi−loops.
#
#   −−noPS  Do not produce postscript drawing of the mfe structure. (default=off)
#
# We use -p to produce the centroid structure (better prediction usually:
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1370799/) . We use -p0 to skip
# the calculation of the pair probabilities, saves time. We use -d2 for -p.

cmd = 'RNAfold -p0 -d2 --noLP --noPS --infile=./Data/Structures/references.fasta --outfile=RNAfold_ref_structures.ifold'
subprocess.run(cmd, shell = True)

# moving the produced .ifold file to the Structures directory

original = './RNAfold_ref_structures.ifold'
target = './Data/Structures/RNAfold_ref_structures.ifold'
shutil.move(original, target)
        
# extracting the centroid structure from the .ifold file
    
with open('./Data/Structures/RNAfold_ref_structures.ifold', 'r') as file:
    line_list = file.readlines()
    for count in range(0, len(line_list)):
        line = line_list[count]
        if line.startswith('>'):
            Type = line.strip()[-3]
            Segment = line.strip()[-1]
        if re.search('[AUTGC]+', line):
            structure = line_list[count+1].strip().split()[0]
            references[Type][Segment]['Ungapped_structure'] = structure

with open('./Data/Metadata/references.json', 'w') as file:
    json.dump(references, file, default = str)
    
#%% ungapped RNA structure figure

URLs = {}
for Type in ['A', 'B']:
    URLs[Type] = {}
    for Segment in list('12345678'):
        Sequence = references[Type][Segment]['Ungapped_sequence']
        Structure = references[Type][Segment]['Ungapped_structure']
        URL = 'http://nibiru.tbi.univie.ac.at/forna/forna.html?id=url/name&sequence='+Sequence+'&structure='+Structure
        URLs[Type][Segment] = URL

Type = 'A'
Segment = '7'

webbrowser.open(URLs[Type][Segment])

#%% Making a fasta to display all structures for each type

# global structures
for Type in ['A', 'B']:
    with open('./Data/Structures/structures_'+Type+'.txt', 'w') as file:
        for Segment in list('12345678'):
            Sequence = references[Type][Segment]['Ungapped_sequence']
            Structure = references[Type][Segment]['Ungapped_structure']           
            file.write('>'+Type+'_'+Segment+'\n'+Sequence+'\n'+Structure+'\n')
            
#%% Making colors files, based on log(mutation count)

for Type in ['A', 'B']:
    with open('./Data/Structures/Colors/colors_'+Type+'.txt', 'w') as all_file:
        for Segment in list('12345678'):
            seg_colors = globals()['df_'+Type+'_'+Segment].sort_values('Position')['Mutation_count'].copy()
            log_colors = np.log10(seg_colors + 1)/max(np.log10(seg_colors + 1))
            all_file.write('>'+Type+'_'+Segment+'\n')
            with open('./Data/Structures/Colors/colors_'+Type+'_'+Segment+'.txt', 'w') as file:
                for color in log_colors:
                    file.write(str(round(color, 2))+'\n')
                    all_file.write(str(round(color, 2))+'\n')

#%% Making histograms

import plotly.express as px
import plotly.io as pio

# setting the renderer to browser
pio.renderers.default = 'browser'
pio.templates.default = "plotly_white"

Type = 'A'
    
Variability = globals()['df_'+Type]['Mutation_count']
# colors = Variability.copy()
# stand_colors = colors/max(colors)


fig = px.histogram(np.log10(Variability + 1)/max(np.log10(Variability + 1)), nbins = 300)
fig.update_layout(
    title_text = 'Variability histogram for Influenza '+Type,
    yaxis = {'title': 'Count'},
    xaxis = {'title': 'Variability'}
    )
fig.show()
        
