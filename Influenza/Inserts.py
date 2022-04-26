#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 25 13:13:23 2021

@author: mennovandamme
"""

### Checking the importance of inserts in the alignment

#%% Dependencies

import json, time, webbrowser, re, subprocess, shutil
# conda install -c anaconda pandas
import pandas as pd
# conda install -c plotly plotly
import plotly.express as px
import plotly.io as pio

# setting the renderer to browser
pio.renderers.default = 'browser'
pio.templates.default = "plotly_white"

#%% Loading the necessary data

with open('./Data/Metadata/metadata_human.json') as file:
    metadata_human = json.load(file)
    
with open('./Data/Metadata/references.json') as file:
    references = json.load(file)

for Type in ['A', 'B']:
    with open('./Data/Dataframes/Complete/complete_df_'+Type+'.csv') as file:
        globals()['complete_df_'+Type] = pd.read_csv(file)
    for Segment in list('12345678'):
        with open('./Data/Dataframes/Complete/complete_df_'+Type+'_'+Segment+'.csv') as file:
            globals()['complete_df_'+Type+'_'+Segment] = pd.read_csv(file)
        
for Type in ['A', 'B']:
    with open('./Data/Dataframes/Ungapped/df_'+Type+'.csv') as file:
        globals()['df_'+Type] = pd.read_csv(file)
    for Segment in list('12345678'):
        with open('./Data/Dataframes/Ungapped/df_'+Type+'_'+Segment+'.csv') as file:
            globals()['df_'+Type+'_'+Segment] = pd.read_csv(file)
            
#%% Gap percentage per type bar chart

for Type in ['A', 'B']:
    with open('./Data/Dataframes/Complete/complete_df_'+Type+'.csv') as file:
        globals()['complete_df_'+Type] = pd.read_csv(file)

Type = 'A'

Color = ['#000000']*len(globals()['complete_df_'+Type]['Position'])

# making animated Graph
fig = px.bar(data_frame = globals()['complete_df_'+Type], 
             x = 'Position', y = 'Gap_percentage', 
             range_y = [0, 105],
             color = Color,
             color_discrete_map = "identity")

fig.update_traces(marker_color = 'black', 
                  marker_line_color = 'black',
                  marker_line_width = 0.05)

# changing axis titles
fig.update_layout(
    title_text = 'Influenza '+Type,
    yaxis = {'title': 'Gap percentage'},
    xaxis = {'title': 'Position'})

# showing the Graph in browser
fig.show()    

#%% Calculating total sequence lengths

for Type in ['A', 'B']:
    tot_length = 0
    for Segment in list('12345678'):
        length = len(references[Type][Segment]['Sequence'])
        tot_length += length
    print(Type+'\t'+str(tot_length))
        
#%% Histogram of gap percentage per position

for Type in ['A', 'B']:
    with open('./Data/Dataframes/Complete/complete_df_'+Type+'.csv') as file:
        globals()['complete_df_'+Type] = pd.read_csv(file)

Type = 'A'
    
gap_percentages = globals()['complete_df_'+Type]['Gap_percentage']

fig = px.histogram(gap_percentages, nbins = 100)
fig.update_layout(
    title_text = 'Gap percentage histogram for Influenza '+Type,
    yaxis = {'title': 'Position count'},
    xaxis = {'title': 'Gap percentage'}
    )

fig.show()

#%% Calculate amount of sequences with rare inserts per segment

for Type in ['A', 'B']:
    for Segment in list('12345678'):
        with open('./Data/Dataframes/Complete/complete_df_'+Type+'_'+Segment+'.csv') as file:
            globals()['complete_df_'+Type+'_'+Segment] = pd.read_csv(file)
                
def insert_counter(Type, Segment):
    total = 0
    count = 0
    complete_df = globals()['complete_df_'+Type+'_'+Segment]
    gap_positions = []
    for i in range(0, len(complete_df)):
        gap_percentage = complete_df.loc[i, 'Gap_percentage']
        if gap_percentage > 80:
            position = complete_df.loc[i, 'Position']
            gap_positions.append(position)
    for ID in metadata_human[Type][Segment]:
        total += 1
        seq = metadata_human[Type][Segment][ID]['Aligned_sequence']
        for pos in gap_positions:
            if seq[pos-1] != '-':
                count += 1
                break
    print(Type+'_'+Segment+': '+str(count)+' out of '+str(total)+' have inserts')    


t2 = time.time()
for Type in ['A', 'B']:
    for Segment in list('12345678'):
        insert_counter(Type, Segment)
t3 = time.time()
print('Calculating the amount of sequences with rare inserts took:', round((t3-t2), 1), 'seconds')  

#%% Calculate amount of sequences with rare inserts per type

def insert_counter2(Type):
    total = 0
    count = 0
    for Segment in list('12345678'):
        complete_df = globals()['complete_df_'+Type+'_'+Segment]
        gap_positions = []
        for i in range(0, len(complete_df)):
            gap_percentage = complete_df.loc[i, 'Gap_percentage']
            if gap_percentage > 90:
                position = complete_df.loc[i, 'Position']
                gap_positions.append(position)
        for ID in metadata_human[Type][Segment]:
            total += 1
            seq = metadata_human[Type][Segment][ID]['Aligned_sequence']
            for pos in gap_positions:
                if seq[pos-1] != '-':
                    count += 1
                    break
    print(Type+': '+str(count)+' out of '+str(total)+' have inserts')    


t2 = time.time()
for Type in ['A', 'B']:
        insert_counter2(Type)
t3 = time.time()
print('Calculating the amount of sequences with rare inserts took:', round((t3-t2), 1), 'seconds')                

