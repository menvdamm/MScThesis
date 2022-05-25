#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 25 14:02:03 2022

@author: mennovandamme
"""

#### Genome viewer for conservation graphs

#%% Dependencies

import os, subprocess
from Bio import SeqIO
import pandas as pd
import numpy as np
from plotly.subplots import make_subplots
import plotly.graph_objects as go
import plotly.io as pio

# setting the renderer to browser
pio.renderers.default = 'browser'
pio.templates.default = "plotly_white"

#%% Directories

if not os.path.isdir('./Data/Genome'): os.mkdir('./Data/Genome')

#%% Alignment

CDS_file = './Data/Genome/CDS.fasta'
consensus_file = './Data/Structure/SARSCoV2_consensus.fasta'
output_file = './Data/Genome/SARSCoV2.fasta'
cmd = 'mafft --auto  --preservecase --keeplength --addfragments '+CDS_file+' --reorder --thread -1 '+consensus_file+' > '+output_file
subprocess.run(cmd, shell = True)

#%% Extract positions of coding sequences

begin_pos = []
end_pos = []
for seq_record in SeqIO.parse('./Data/Genome/SARSCoV2.fasta', 'fasta'):
    ID = seq_record.id
    if ID != 'SARS-CoV-2':
        seq = str(seq_record.seq)
        first = True
        for i in range(len(seq)):
            if seq[i] != '-' and first == True:
                begin_pos.append(i+1)
                first = False
            if first == False and seq[i] == '-':
                end_pos.append(i)
                break
                
CDS_df = pd.DataFrame({'Begin_position': begin_pos,
                       'End_position': end_pos,
                       'Protein': ['ORF1ab', 'ORF1a', 'S', 'ORF3a', 'E', 'M', 'ORF6', 'ORF7a', 'ORF7b', 'ORF8', 'N', 'ORF10']})            

CDS_df.to_csv('./Data/Genome/CDS_df.csv', index = False)

#%% Loading the necessary data

with open('./Data/Dataframes/df.csv') as file:
    df = pd.read_csv(file)

with open('./Data/Genome/CDS_df.csv') as file:
    CDS_df = pd.read_csv(file)
    
#%% Make figure

metric = 'Shannon_entropy'

fig = make_subplots(rows = 2, 
                    cols = 1, 
                    shared_xaxes = True, 
                    row_heights = [0.5, 0.5], 
                    vertical_spacing = 0.01,
                    )

fig.add_trace(
    go.Bar(x = df['Position'], 
           y = df[metric],
           marker_color = 'black',
           marker_line_color = 'black',
           marker_line_width = 0.5,
           showlegend = False,
           name = ' '.join(metric.split('_')),
           ),
    row=1, col=1
)

fig.add_trace(
    go.Bar(base = CDS_df['Begin_position'],
           x = CDS_df['End_position'].to_numpy() - CDS_df['Begin_position'].to_numpy() + 1, 
           y = [1, 2, 2, 1, 2, 2, 1, 1, 1, 1, 2, 1], 
           name = '', 
           showlegend = False,
           orientation = 'h',
           customdata = CDS_df['Protein'],
           texttemplate = "%{customdata}",
           textposition = "inside",
           textangle = 0,
           textfont_color = "black",
           ),
    row=2, col=1
)

fig.update_layout(barmode='stack')

fig.update_xaxes(title_text = 'Position', row = 2)
fig.update_yaxes(title_text = ' '.join(metric.split('_')), row = 1, range = [0, max(df[metric])])
fig.update_yaxes(visible = False, showticklabels = False, row = 2)             

fig.show()







