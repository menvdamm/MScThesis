#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 24 18:07:53 2022

@author: mennovandamme
"""

### RNA secondary structure

#%% Dependencies

import subprocess, os, shutil
import pandas as pd
import numpy as np

#%% Directories

if not os.path.isdir('./Data/Structure'): os.mkdir('./Data/Structure')
if not os.path.isdir('./Data/Structure/Consensus'): os.mkdir('./Data/Structure/Consensus')
if not os.path.isdir('./Data/Structure/Dotbracket'): os.mkdir('./Data/Structure/Dotbracket')
if not os.path.isdir('./Data/Structure/Colors'): os.mkdir('./Data/Structure/Colors')
if not os.path.isdir('./Data/Structure/Element'): os.mkdir('./Data/Structure/Element')

#%% Loading the necessary data

for Type in ['A', 'B']:
    with open('./Data/Dataframes/df_'+Type+'.csv') as file:
        globals()['df_'+Type] = pd.read_csv(file)
    for Segment in list('12345678'): 
        with open('./Data/Dataframes/Ungapped/df_'+Type+'_'+Segment+'.csv') as file:
            globals()['df_'+Type+'_'+Segment] = pd.read_csv(file)

#%% Making consensus sequence 

# per segment
for Type in ['A', 'B']:
    for Segment in list('12345678'): 
        with open('./Data/Structure/Consensus/'+Type+'_'+Segment+'.fasta', 'w') as f:
            df = globals()['df_'+Type+'_'+Segment]
            seq = ''
            for row in df.sort_values('Position').iterrows():
                percentages = {'A': row[1]['A_percentage'],
                               'U': row[1]['T_percentage'],
                               'G': row[1]['G_percentage'],
                               'C': row[1]['C_percentage']}
                most_occuring = max(percentages, key = percentages.get)
                seq += most_occuring
            f.write('>'+Type+'_'+Segment+'\n'+seq)

# per type
for Type in ['A', 'B']:
    sequence = ''
    for Segment in list('12345678'): 
        with open('./Data/Structure/Consensus/'+Type+'_'+Segment+'.fasta', 'r') as f:
            lines = f.readlines()
            seq = lines[1].strip()
            sequence += seq
    with open('./Data/Structure/Consensus/'+Type+'.fasta', 'w') as f:
        f.write('>'+Type+'\n'+sequence)
            
#%% RNA secondary structure prediction

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
 
for Type in ['A', 'B']:
    for Segment in list('12345678'):   
        input_file = './Data/Structure/Consensus/'+Type+'_'+Segment+'.fasta'
        output_file = Type+'_'+Segment+'.ifold'
        cmd = 'RNAfold -p0 -d2 --noLP --noPS --infile='+input_file+' --outfile='+output_file
        subprocess.run(cmd, shell = True)
        # moving the produced .ifold file to the Dotbracket directory
        original = Type+'_'+Segment+'.ifold'
        target = './Data/Structure/Dotbracket/'+original
        shutil.move(original, target)
        # extracting the centroid structure from the .ifold file
        with open('./Data/Structure/Dotbracket/'+Type+'_'+Segment+'.ifold', 'r') as f:
            lines = f.readlines()
            sequence = lines[1].strip()
            structure = lines[2].strip().split(' ')[0]
        # making a file with sequence and structure
        with open('./Data/Structure/Dotbracket/'+Type+'_'+Segment+'.txt', 'w') as file:         
            file.write('>'+Type+'_'+Segment+'\n'+sequence+'\n'+structure)

# combine into one file per type
for Type in ['A', 'B']:
    with open('./Data/Structure/Dotbracket/'+Type+'.txt', 'w') as type_file:
        with open('./Data/Structure/Dotbracket/'+Type+'_'+Segment+'.txt', 'r') as segment_file:
            lines = segment_file.readlines()
            header = lines[0].strip()
            sequence = lines[1].strip()
            structure = lines[2].strip()          
            type_file.write(header+'\n'+sequence+'\n'+structure+'\n')       

#%% Making colors files for visualization in forna, based on log(mutability)

for Type in ['A', 'B']:
    with open('./Data/Structure/Colors/'+Type+'.txt', 'w') as type_file:
        for Segment in list('12345678'):
            mutability = globals()['df_'+Type+'_'+Segment].sort_values('Position')['Mutability'].copy()
            colors = np.log10(mutability*1000 + 1)/max(np.log10(mutability*10000 + 1))
            type_file.write('>'+Type+'_'+Segment+'\n')
            with open('./Data/Structure/Colors/'+Type+'_'+Segment+'.txt', 'w') as segment_file:
                for color in colors:
                    segment_file.write(str(round(color, 2))+'\n')
                    type_file.write(str(round(color, 2))+'\n')

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

for Type in ['A', 'B']:
    for Segment in list('12345678'):   
        input_file = './Data/Structure/Dotbracket/'+Type+'_'+Segment+'.txt'
        output_file = './Data/Structure/Element/'+Type+'_'+Segment+'_'
        cmd = 'rnaConvert.py '+input_file+' -T element_string --filename '+output_file
        subprocess.run(cmd, shell = True)

#%% Adding structure element to conservation dataframes

for Type in ['A', 'B']:
    elements = ''
    for Segment in list('12345678'):
        with open('./Data/Structure/Element/'+Type+'_'+Segment+'_001.element_string') as file:
             ele = file.readlines()[1].strip()
             elements += ele
             globals()['df_'+Type+'_'+Segment] = globals()['df_'+Type+'_'+Segment].sort_values('Position')
             globals()['df_'+Type+'_'+Segment]['Element'] = list(ele)
             globals()['df_'+Type+'_'+Segment]['Element'] = globals()['df_'+Type+'_'+Segment]['Element'].replace({"f": "5' unpaired",                
                                                                                                                "i": "interior loop & bulge",
                                                                                                                "s": "stem",
                                                                                                                "m": "multiloop",
                                                                                                                "h": "hairpin loop",
                                                                                                                "t": "3' unpaired"
                                                                                                                })
             globals()['df_'+Type+'_'+Segment].to_csv('./Data/Dataframes/Ungapped/df_'+Type+'_'+Segment+'.csv', index = False)
    globals()['df_'+Type] = globals()['df_'+Type].sort_values('Position')
    globals()['df_'+Type]['Element'] = list(elements)
    globals()['df_'+Type]['Element'] = globals()['df_'+Type]['Element'].replace({"f": "5' unpaired",                
                                                                                "i": "interior loop & bulge",
                                                                                "s": "stem",
                                                                                "m": "multiloop",
                                                                                "h": "hairpin loop",
                                                                                "t": "3' unpaired"
                                                                                })
    globals()['df_'+Type].to_csv('./Data/Dataframes/Ungapped/df_'+Type+'.csv', index = False)

#%% Statistical testing

## chi-squared goodness-of-fit test versus a random distribution

# for Type in ['A', 'B']:
#     # generating random data
#     random_mutation_count = globals()['df_'+Type]['Mutation_count'].sample(frac=1).reset_index(drop=True)
#     iterations = 5000
#     count = 1
#     while count < iterations:
#         random_mutation_count += globals()['df_'+Type]['Mutation_count'].sample(frac=1).reset_index(drop=True)
#         count += 1
#     random_data = pd.DataFrame({'Random_mutation_count': random_mutation_count,
#                                 'Element': globals()['df_'+Type]['Element']})   
#     df = pd.DataFrame()
#     df['Element_count'] = pd.DataFrame(globals()['df_'+Type][{'Element', 'Position'}].groupby(['Element']).count())
#     # Mutation_count * iterations to get the same total count in Mutation_count & Random_mutation_count
#     df['Mutation_count'] = pd.DataFrame(globals()['df_'+Type][{'Element', 'Mutation_count'}].groupby(['Element']).sum()) * iterations
#     df['Random_mutation_count'] = pd.DataFrame(random_data.groupby(['Element']).sum())
#     globals()['stat_df_'+Type] = df.copy()
    
# chisquare(f_obs = stat_df_A['Mutation_count'], f_exp = stat_df_A['Random_mutation_count'])

# chisquare(f_obs = stat_df_B['Mutation_count'], f_exp = stat_df_B['Random_mutation_count'])
     