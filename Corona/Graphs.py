#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 12 13:11:11 2021

@author: mennovandamme
"""

#### Displaying the data

#%% Dependencies

import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.io as pio


# setting the renderer to browser
pio.renderers.default = 'browser'
pio.templates.default = "plotly_white"

#%% Loading the necessary data

with open('./Data/Dataframes/complete_df.csv') as file:
    complete_df = pd.read_csv(file)
    
with open('./Data/Dataframes/df.csv') as file:
    df = pd.read_csv(file)
    
with open('./Data/Dataframes/score_df.csv') as file:
    score_df = pd.read_csv(file)

with open('./Data/Genome/CDS_df.csv') as file:
    CDS_df = pd.read_csv(file)

#%% Gap percentage

Color = ['#000000']*len(complete_df)

fig = px.bar(data_frame = complete_df, 
             x = 'Position', y = 'Gap_percentage',
             range_y = [0, 100],
             color = Color,
             color_discrete_map = "identity"
             )

fig.update_traces(marker_color = 'black', 
                  marker_line_color = 'black',
                  marker_line_width = 0.05)

fig.update_layout(
#    title_text = 'Gap percentage per position for the SARS-CoV-2 alignment',
    yaxis = {'title': 'Mutability'},
    xaxis = {'title': 'Position'}
    )

fig.show()

#%% Bar charts

Color = ['#000000']*len(df)


fig = px.bar(data_frame = df, 
             x = 'Position', y = 'Mutability',
             range_y = [0, max(df['Mutability'])],
             color = Color,
             color_discrete_map = "identity"
             )

fig.update_traces(marker_color = 'black', 
                  marker_line_color = 'black',
                  marker_line_width = 0.05)


fig.update_layout(
#    title_text = 'Mutability per position for SARS-CoV-2',
    yaxis = {'title': 'Mutability'},
    xaxis = {'title': 'Position'}
    )

fig.show()


fig = px.bar(data_frame = df, 
             x = 'Position', y = 'Shannon_entropy',
             range_y = [0, max(df['Shannon_entropy'])],
             color = Color,
             color_discrete_map = "identity"
             )

fig.update_traces(marker_color = 'black', 
                  marker_line_color = 'black',
                  marker_line_width = 0.05)

fig.update_layout(
#    title_text = 'Shannon entropy per position for SARS-CoV-2',
    yaxis = {'title': 'Shannon entropy'},
    xaxis = {'title': 'Position'}
    )

fig.show()

#%% Histograms
    
fig = px.histogram(df['Mutability'], nbins = 1000)

fig.update_layout(
#    title_text = 'Mutability histogram for SARS-CoV-2',
    yaxis = {'title': 'Number of positions'},
    xaxis = {'title': 'Mutability'},
    showlegend=False
    )

fig.show()
    

fig = px.histogram(df['Shannon_entropy'], nbins = 1000)

fig.update_layout(
#    title_text = 'Shannon entropy histogram for SARS-CoV-2',
    yaxis = {'title': 'Number of positions'},
    xaxis = {'title': 'Shannon entropy'},
    showlegend=False
    )

fig.show()

#%% Correlation between scores

fig = px.scatter(df, x="Mutability", y="Shannon_entropy") #, trendline="ols")

fig.update_layout(
#    title_text = 'Shannon entropy in function of Mutability per position for SARS-CoV-2',
    yaxis = {'title': 'Shannon entropy'},
    xaxis = {'title': 'Mutability'}
    )

fig.show()

#%% Bar chart with most conserved 30b regions

fig = go.Figure()

fig.add_trace(
    go.Bar(
        x = df['Position'],
        y = df['Shannon_entropy'],
        marker_color = 'black',
        marker_line_color = 'black',
        marker_line_width = 0.05,
        showlegend = False))
        
sdf = score_df.sort_values('Shannon_entropy')[0:10]
for index, row in sdf.iterrows():
    fig.add_trace(
        go.Scatter(
            x = [row['Begin_position'], row['End_position']],
            y = [-0.1, -0.1],
            mode = 'lines',
#                name = '30 bp',
            line = dict(color = 'blue', width = 6),
            showlegend = False))

fig.update_layout(
#    title_text = 'Shannon entropy per position for SARS-CoV-2',
    yaxis = {'title': 'Shannon entropy'},
    xaxis = {'title': 'Position'}
    )

fig.show()


fig = go.Figure()

fig.add_trace(
    go.Bar(
        x = df['Position'],
        y = df['Shannon_entropy'],
        marker_color = 'black',
        marker_line_color = 'black',
        marker_line_width = 0.05,
        showlegend = False))
        

sdf = score_df.sort_values('Mutability')[0:10]
for index, row in sdf.iterrows():
    fig.add_trace(
        go.Scatter(
            x = [row['Begin_position'], row['End_position']],
            y = [-0.1, -0.1],
            mode = 'lines',
#           name = '30 bp',
            line = dict(color = 'blue', width = 6),
            showlegend = False))

fig.update_layout(
#    title_text = 'Mutability per position for SARS-CoV-2',
    yaxis = {'title': 'Mutability'},
    xaxis = {'title': 'Position'}
    )

fig.show()

#%% Bar charts per 30b

Color = ['#000000']*len(score_df)

fig = px.bar(data_frame = score_df, 
             x = 'Begin_position', y = 'Mutability',
             range_y = [0, max(score_df['Mutability'])],
             color = Color,
             color_discrete_map = "identity"
             )

fig.update_traces(marker_color = 'black', 
                  marker_line_color = 'black',
                  marker_line_width = 0.05)

fig.update_layout(
#    title_text = 'Mutability per 30b for SARS-CoV-2',
    yaxis = {'title': 'Mutability'},
    xaxis = {'title': 'Begin position'}
    )

fig.show()


fig = px.bar(data_frame = score_df, 
             x = 'Begin_position', y = 'Shannon_entropy',
             range_y = [0, max(score_df['Shannon_entropy'])],
             color = Color,
             color_discrete_map = "identity"
             )

fig.update_traces(marker_color = 'black', 
                  marker_line_color = 'black',
                  marker_line_width = 0.05)

fig.update_layout(
#    title_text = 'Shannon entropy per 30b for SARS-CoV-2',
    yaxis = {'title': 'Shannon entropy'},
    xaxis = {'title': 'Begin position'}
    )

fig.show()

#%% Histograms per 30 bp
    
fig = px.histogram(score_df['Mutability'], nbins = 1000)

fig.update_layout(
#    title_text = 'Mutability histogram per 30b for SARS-CoV-2',
    yaxis = {'title': 'Number of candidates'},
    xaxis = {'title': 'Mutability'},
    showlegend=False
    )

fig.show()
    

fig = px.histogram(score_df['Shannon_entropy'], nbins = 1000)

fig.update_layout(
#    title_text = 'Shannon entropy histogram per 30b for SARS-CoV-2',
    yaxis = {'title': 'Number of candidates'},
    xaxis = {'title': 'Shannon entropy'},
    showlegend=False
    )

fig.show()

#%% Correlation between scores per 30 bp

fig = px.scatter(score_df, x="Mutability", y="Shannon_entropy") #, trendline="ols")

fig.update_layout(
#    title_text = 'Shannon entropy in function of Mutability per 30b for SARS-CoV-2',
    yaxis = {'title': 'Shannon entropy'},
    xaxis = {'title': 'Mutability'}
    )

fig.show()

#%% RNA structure color histogram


colors = df.sort_values('Position')['Mutability'].copy()
log_colors = np.log10(colors*1000 + 1)/max(np.log10(colors*10000 + 1))

fig = px.histogram(log_colors, nbins = 300)

fig.update_layout(
#    title_text = 'Color histogram for SARS-CoV-2',
    yaxis = {'title': 'Count'},
    xaxis = {'title': 'Color'}
    )

fig.show()

#%% Element structure Box plot

fig = px.box(df, x="Element", y="Shannon_entropy", category_orders={'Element': ["5' unpaired", "stem", "hairpin loop", "interior loop & bulge", "multiloop", "3' unpaired"]})

fig.update_layout(
#    title_text = 'RNA structure element Shannon entropy boxplot for for SARS-CoV-2',
    yaxis = {'title': 'Shannon entropy'},
    xaxis = {'title': 'Structure element'}
    )

fig.show()

#%% Conservation bar chart with genome viewer

# choose what you want to display

amount_of_sequences = 10

metric = 'Mutability' #'Shannon_entropy'

if metric == 'Shannon_entropy':
    y_cons = [2, 1, 1, 1, 1, 2, 2, 2, 1, 2] + [1]*1000
elif metric == 'Mutability':
    y_cons = [2, 1, 1, 1, 1, 1, 2, 2, 1, 1] + [1]*1000

windowsize = 1 #30

if windowsize == 1:
    x = df['Position']
    y = df[metric]
    x_axis = 'Position'
    line_thickness = 1
elif windowsize == 30:
    x = score_df['Begin_position']
    y = score_df[metric]
    x_axis = 'Begin position'
    line_thickness = 0.5

fig = make_subplots(rows = 2, 
                    cols = 1, 
                    shared_xaxes = True, 
                    row_heights = [0.90, 0.1], 
                    vertical_spacing = 0.01,
                    )

fig.add_trace(
    go.Bar(x = x, 
           y = y,
           marker_color = 'black',
           marker_line_color = 'black',
           marker_line_width = line_thickness,
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
           textposition = 'inside',
           insidetextanchor = 'middle',
           textangle = 0,
           ),
    row=2, col=1
)

small_score_df = score_df.sort_values(metric)[0:amount_of_sequences]

fig.add_trace(
    go.Bar(base = small_score_df['Begin_position'],
           x = small_score_df['End_position'].to_numpy() - small_score_df['Begin_position'].to_numpy() + 1, 
           y = y_cons, 
           name = '', 
           showlegend = False,
           orientation = 'h',
           marker_color = 'blue',
           marker_line_color = 'blue',
           ),
    row=2, col=1
)

fig.update_layout(barmode='stack', uniformtext_minsize=6, uniformtext_mode='show')

fig.update_xaxes(title_text = x_axis, row = 2)
fig.update_yaxes(title_text = ' '.join(metric.split('_')), row = 1, range = [0, max(df[metric])], fixedrange = True)
fig.update_yaxes(visible = False, showticklabels = False, row = 2, fixedrange = True)             

fig.show()