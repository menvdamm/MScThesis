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

for Type in ['A', 'B']:
    with open('./Data/Dataframes/df_'+Type+'.csv') as file:
        globals()['df_'+Type] = pd.read_csv(file)
    with open('./Data/Dataframes/score_df_'+Type+'.csv') as file:
        globals()['score_df_'+Type] = pd.read_csv(file)
    with open('./Data/Genome/CDS_df_'+Type+'.csv') as file:
        globals()['CDS_df_'+Type] = pd.read_csv(file, keep_default_na=False)
    for Segment in list('12345678'): 
        with open('./Data/Dataframes/Ungapped/df_'+Type+'_'+Segment+'.csv') as file:
            globals()['df_'+Type+'_'+Segment] = pd.read_csv(file)
        with open('./Data/Dataframes/Complete/df_'+Type+'_'+Segment+'.csv') as file:
            globals()['complete_df_'+Type+'_'+Segment] = pd.read_csv(file)
        

#%% Choose type & segment to display data for

Type = 'A'

Segment = '1'

complete_df = globals()['complete_df_'+Type+'_'+Segment]
# df = globals()['df_'+Type+'_'+Segment]
df = globals()['df_'+Type]
score_df = globals()['score_df_'+Type]
CDS_df = globals()['CDS_df_'+Type]

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

# changing axis titles
fig.update_layout(
#    title_text = 'Gap percentage per position for the Influenza '+Type+' Segment '+Segment+' alignment',
    yaxis = {'title': 'Mutability'},
    xaxis = {'title': 'Position'}
    )

# showing the Graph in browser
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

# changing axis titles
fig.update_layout(
#    title_text = 'Mutability per position for Influenza '+Type+' Segment '+Segment,
    yaxis = {'title': 'Mutability'},
    xaxis = {'title': 'Position'}
    )

# showing the Graph in browser
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

# changing axis titles
fig.update_layout(
#    title_text = 'Shannon entropy per position for Influenza '+Type+' Segment '+Segment,
    yaxis = {'title': 'Shannon entropy'},
    xaxis = {'title': 'Position'}
    )

# showing the Graph in browser
fig.show()

#%% Histograms
    
fig = px.histogram(df['Mutability'], nbins = 1000)

fig.update_layout(
#    title_text = 'Mutability histogram for Influenza '+Type+' Segment '+Segment,
    yaxis = {'title': 'Position count'},
    xaxis = {'title': 'Mutability'}
    )

fig.show()
    

fig = px.histogram(df['Shannon_entropy'], nbins = 1000)

fig.update_layout(
#    title_text = 'Shannon entropy histogram for Influenza '+Type+' Segment '+Segment,
    yaxis = {'title': 'Position count'},
    xaxis = {'title': 'Shannon entropy'}
    )

fig.show()

#%% Correlation between scores

fig = px.scatter(df, x="Mutability", y="Shannon_entropy") #, trendline="ols")

fig.update_layout(
#    title_text = 'Shannon entropy in function of Mutability per position for Influenza '+Type+' Segment '+Segment,
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
#    title_text = 'Shannon entropy per position for Influenza '+Type+' Segment '+Segment,
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
#                name = '30 bp',
            line = dict(color = 'blue', width = 6),
            showlegend = False))

fig.update_layout(
#    title_text = 'Mutability per position for Influenza '+Type+' Segment '+Segment,
    yaxis = {'title': 'Mutability'},
    xaxis = {'title': 'Position'}
    )

fig.show()

#%% Bar charts per 30b

Color = ['#000000']*len(score_df)

# making animated Graph
fig = px.bar(data_frame = score_df, 
             x = 'Begin_position', y = 'Mutability',
             range_y = [0, max(score_df['Mutability'])],
             color = Color,
             color_discrete_map = "identity"
             )

fig.update_traces(marker_color = 'black', 
                  marker_line_color = 'black',
                  marker_line_width = 0.05)

# changing axis titles
fig.update_layout(
    title_text = 'Mutability per 30b for Influenza '+Type+' Segment '+Segment,
    yaxis = {'title': 'Mutability'},
    xaxis = {'title': 'Begin_position'}
    )

# showing the Graph in browser
fig.show()


# making animated Graph
fig = px.bar(data_frame = score_df, 
             x = 'Begin_position', y = 'Shannon_entropy',
             range_y = [0, max(score_df['Shannon_entropy'])],
             color = Color,
             color_discrete_map = "identity"
             )

fig.update_traces(marker_color = 'black', 
                  marker_line_color = 'black',
                  marker_line_width = 0.05)

# changing axis titles
fig.update_layout(
    title_text = 'Shannon entropy per 30b for Influenza '+Type+' Segment '+Segment,
    yaxis = {'title': 'Shannon entropy'},
    xaxis = {'title': 'Begin position'}
    )

# showing the Graph in browser
fig.show()

#%% Histograms per 30 bp
    
fig = px.histogram(score_df['Mutability'], nbins = 300)

fig.update_layout(
    title_text = 'Mutability histogram per 30b for Influenza '+Type+' Segment '+Segment,
    yaxis = {'title': 'Position count'},
    xaxis = {'title': 'Mutability'}
    )

fig.show()
    

fig = px.histogram(score_df['Shannon_entropy'], nbins = 300)

fig.update_layout(
    title_text = 'Shannon entropy histogram per 30b for Influenza '+Type+' Segment '+Segment,
    yaxis = {'title': 'Position count'},
    xaxis = {'title': 'Shannon entropy'}
    )

fig.show()

#%% Correlation between scores per 30 bp

fig = px.scatter(score_df, x="Mutability", y="Shannon_entropy") #, trendline="ols")

fig.update_layout(
    title_text = 'Shannon entropy in function of Mutability per 30b for Influenza '+Type+' Segment '+Segment,
    yaxis = {'title': 'Shannon entropy'},
    xaxis = {'title': 'Mutability'}
    )

fig.show()

#%% RNA structure color histogram


colors = df.sort_values('Position')['Mutability'].copy()
log_colors = np.log10(colors*1000 + 1)/max(np.log10(colors*10000 + 1))

fig = px.histogram(log_colors, nbins = 300)

fig.update_layout(
    title_text = 'Color histogram for Influenza '+Type+' Segment '+Segment,
    yaxis = {'title': 'Count'},
    xaxis = {'title': 'Color'}
    )

fig.show()

#%% Element structure Box plot

fig = px.box(df, x="Element", y="Shannon_entropy")

fig.update_layout(
    title_text = 'RNA structure element Shannon entropy boxplot for for Influenza '+Type+' Segment '+Segment,
    yaxis = {'title': 'Shannon entropy'},
    xaxis = {'title': 'Structure element'}
    )

fig.show()

#%% Conservation bar chart with genome viewer

metric = 'Shannon_entropy'

if Type == 'A':
    y_CDS = [2, 2, 1, 1, 1, 2, 1, 1, 1, 2, 1, 1, 2]
    y_cons = [2]*10
elif Type == 'B':
    y_CDS = [2, 2, 1, 1, 1, 2, 1, 1, 1, 2, 1, 1, 2]
    y_cons = [2]*10

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
           marker_line_width = 0.25,
           showlegend = False,
           name = ' '.join(metric.split('_')),
           ),
    row=1, col=1
)

fig.add_trace(
    go.Bar(base = CDS_df['Begin_position'],
           x = CDS_df['End_position'].to_numpy() - CDS_df['Begin_position'].to_numpy() + 1, 
           y = y_CDS, 
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

small_score_df = score_df.sort_values(metric).copy()[:][0:10]

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

fig.update_layout(barmode='stack', uniformtext_minsize=8, uniformtext_mode='show')

fig.update_xaxes(title_text = 'Position', row = 2)
fig.update_yaxes(title_text = ' '.join(metric.split('_')), row = 1, range = [0, max(df[metric])], fixedrange = True)
fig.update_yaxes(visible = False, showticklabels = False, row = 2, fixedrange = True)             

fig.show()