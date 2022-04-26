#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  2 17:47:36 2021

@author: mennovandamme
"""

### Secondary RNA structure statistical analysis (SSA)

#%% Dependencies

# conda activate Masterproef
# conda install pip
# /opt/anaconda3/envs/Masterproef/bin/pip install forgi
import forgi
import json, subprocess, os, time
import pandas as pd
import numpy as np
import statsmodels
from scipy.stats import chisquare
import plotly.express as px

#%% Loading data

with open('./Data/Metadata/references.json') as file:
    references = json.load(file)
    
for Type in ['A', 'B']:
    with open('./Data/Dataframes/Ungapped/df_'+Type+'.csv') as file:
        globals()['df_'+Type] = pd.read_csv(file)

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
    cmd = 'rnaConvert.py ./Data/Structures/structures_'+Type+'.txt -T element_string --filename ./Data/Structures/element_structure_'
    subprocess.run(cmd, shell = True)
    for Segment in list('12345678'):
        old_file = os.path.join('./Data/Structures/', 'element_structure_00'+Segment+'.element_string')
        new_file = os.path.join('./Data/Structures/', 'element_structure_'+Type+'_'+Segment+'.txt')
        os.rename(old_file, new_file)

#%% Adding element to dataframe

for Type in ['A', 'B']:
    elements = ''
    for Segment in list('12345678'):
        with open('./Data/Structures/element_structure_'+Type+'_'+Segment+'.txt') as file:
            references[Type][Segment]['Element_structure'] = file.readlines()[1].strip()
            elements += references[Type][Segment]['Element_structure']
    globals()['df_'+Type] = globals()['df_'+Type].sort_values('Position')
    globals()['df_'+Type]['Element'] = list(elements)
    globals()['df_'+Type].to_csv('./Data/Dataframes/Ungapped/df_'+Type+'.csv', index = False)
                
#%% Box plot

Type = 'A'

df = globals()['df_'+Type].copy()

df['Element'] = df['Element'].replace({"f": "5' unpaired",                
          "i": "interior loop & bulge",
          "s": "stem",
          "m": "multiloop",
          "h": "hairpin loop",
          "t": "3' unpaired"
          })

#fig = px.box(globals()['df_'+Type], x="Element", y="Variability")
#fig.show()

#fig = px.box(globals()['df_'+Type], x="Element", y="Mutation_count")
#fig.show()

fig = px.box(df, x="Element", y="Shannon_entropy")

fig.update_layout(
    title_text = 'RNA structure element Shannon entropy boxplot for Influenza '+Type,
    yaxis = {'title': 'Shannon entropy'},
    xaxis = {'title': 'Structure element'}
    )

fig.show()

#%% Statistical testing

## chi-squared goodness-of-fit test versus a random distribution

for Type in ['A', 'B']:
    # generating random data
    random_mutation_count = globals()['df_'+Type]['Mutation_count'].sample(frac=1).reset_index(drop=True)
    iterations = 5000
    count = 1
    while count < iterations:
        random_mutation_count += globals()['df_'+Type]['Mutation_count'].sample(frac=1).reset_index(drop=True)
        count += 1
    random_data = pd.DataFrame({'Random_mutation_count': random_mutation_count,
                                'Element': globals()['df_'+Type]['Element']})   
    df = pd.DataFrame()
    df['Element_count'] = pd.DataFrame(globals()['df_'+Type][{'Element', 'Position'}].groupby(['Element']).count())
    # Mutation_count * iterations to get the same total count in Mutation_count & Random_mutation_count
    df['Mutation_count'] = pd.DataFrame(globals()['df_'+Type][{'Element', 'Mutation_count'}].groupby(['Element']).sum()) * iterations
    df['Random_mutation_count'] = pd.DataFrame(random_data.groupby(['Element']).sum())
    globals()['stat_df_'+Type] = df.copy()
    
chisquare(f_obs = stat_df_A['Mutation_count'], f_exp = stat_df_A['Random_mutation_count'])

chisquare(f_obs = stat_df_B['Mutation_count'], f_exp = stat_df_B['Random_mutation_count'])


