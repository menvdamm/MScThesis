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
from plotly.subplots import make_subplots
import plotly.graph_objects as go
import plotly.io as pio

# setting the renderer to browser
pio.renderers.default = 'browser'
pio.templates.default = "plotly_white"

#%% Loading the necessary data

with open('./Data/Dataframes/df.csv') as file:
    df = pd.read_csv(file)

with open('./Data/Dataframes/score_df.csv') as file:
    score_df = pd.read_csv(file)

with open('./Data/Genome/CDS_df.csv') as file:
    CDS_df = pd.read_csv(file)
    
#%% Make figure

metric = 'Shannon_entropy'

fig = make_subplots(rows = 2, 
                    cols = 1, 
                    shared_xaxes = True, 
                    row_heights = [0.90, 0.1], 
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
           y = [1, 2, 2, 1, 2, 2, 1, 0, 1, 0, 2, 1], 
           name = '', 
           showlegend = False,
           marker_color = 'orange',
           marker_line_color = 'orange',
           opacity = 0.8,
           orientation = 'h',
           customdata = CDS_df['Protein'],
           texttemplate = "%{customdata}",
           textposition = 'outside',
 #          insidetextanchor = 'middle',
           textangle = 0,
           ),
    row=2, col=1
)

small_score_df = score_df.sort_values(metric).copy()[:][0:10]

fig.add_trace(
    go.Bar(base = small_score_df['Begin_position'],
           x = small_score_df['End_position'].to_numpy() - small_score_df['Begin_position'].to_numpy() + 1, 
           y = [2, 1, 1, 1, 1, 2, 2, 2, 1, 2], 
           name = '', 
           showlegend = False,
           orientation = 'h',
           marker_color = 'blue',
           marker_line_color = 'blue',
           ),
    row=2, col=1
)

fig.update_layout(barmode='stack', uniformtext_minsize=8, uniformtext_mode='show')

fig.update_xaxes(title_text = 'Position', row = 2)
fig.update_yaxes(title_text = ' '.join(metric.split('_')), row = 1, range = [0, max(df[metric])], fixedrange = True)
fig.update_yaxes(visible = False, showticklabels = False, row = 2, fixedrange = True)             

fig.show()







