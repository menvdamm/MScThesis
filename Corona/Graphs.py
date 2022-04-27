#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 12 13:11:11 2021

@author: mennovandamme
"""

#### Displaying the animated mutation figures

#%% Dependencies

import pandas as pd
# conda install -c plotly plotly
import plotly.express as px
import plotly.graph_objects as go
# conda install statsmodels
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

#%% Gap percentage

Color = ['#000000']*len(complete_df)

# making animated Graph
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
    title_text = 'Gap percentage per position for the SARS-CoV-2 alignment',
    yaxis = {'title': 'Mutability'},
    xaxis = {'title': 'Position'}
    )

# showing the Graph in browser
fig.show()

#%% Bar charts

Color = ['#000000']*len(df)

# making animated Graph
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
    title_text = 'Mutability per position for SARS-CoV-2',
    yaxis = {'title': 'Mutability'},
    xaxis = {'title': 'Position'}
    )

# showing the Graph in browser
fig.show()


# making animated Graph
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
    title_text = 'Shannon entropy per position for SARS-CoV-2',
    yaxis = {'title': 'Shannon entropy'},
    xaxis = {'title': 'Position'}
    )

# showing the Graph in browser
fig.show()

#%% Histograms
    
fig = px.histogram(df['Mutability'], nbins = 300)

fig.update_layout(
    title_text = 'Mutability histogram for SARS-CoV-2',
    yaxis = {'title': 'Position count'},
    xaxis = {'title': 'Mutability'}
    )

fig.show()
    

fig = px.histogram(df['Shannon_entropy'], nbins = 300)

fig.update_layout(
    title_text = 'Shannon entropy histogram for SARS-CoV-2',
    yaxis = {'title': 'Position count'},
    xaxis = {'title': 'Shannon entropy'}
    )

fig.show()

#%% Correlation between scores

fig = px.scatter(df, x="Mutability", y="Shannon_entropy") #, trendline="ols")

fig.update_layout(
    title_text = 'Shannon entropy in function of Mutability per position for SARS-CoV-2',
    yaxis = {'title': 'Shannon entropy'},
    xaxis = {'title': 'Mutability'}
    )

fig.show()

#%% Histogram with most conserved 39b regions

fig = go.Figure()

fig.add_trace(
    go.Bar(
        x = df['Position'],
        y = df['Shannon_entropy'],
        marker_color = 'black',
        marker_line_color = 'black',
        marker_line_width = 0.05,
        showlegend = False))
        
for var in ['Mutability', 'Shannon_entropy']:
    sdf = score_df.sort_values(var)[0:10]
    for index, row in sdf.iterrows():
        fig.add_trace(
            go.Scatter(
                x = [row['Begin_position'], row['End_position']],
                y = [-0.1, -0.1],
                mode = 'lines',
#                name = '39 bp',
                line = dict(color = 'blue', width = 6),
                showlegend = False))
   
fig.show()

#%% Bar charts per 39b

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
    title_text = 'Mutability per 39b for SARS-CoV-2',
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
    title_text = 'Shannon entropy per 39b for SARS-CoV-2',
    yaxis = {'title': 'Shannon entropy'},
    xaxis = {'title': 'Begin_position'}
    )

# showing the Graph in browser
fig.show()

#%% Histograms
    
fig = px.histogram(score_df['Mutability'], nbins = 300)

fig.update_layout(
    title_text = 'Mutability histogram per 39b for SARS-CoV-2',
    yaxis = {'title': 'Position count'},
    xaxis = {'title': 'Mutability'}
    )

fig.show()
    

fig = px.histogram(score_df['Shannon_entropy'], nbins = 300)

fig.update_layout(
    title_text = 'Shannon entropy histogram per 39b for SARS-CoV-2',
    yaxis = {'title': 'Position count'},
    xaxis = {'title': 'Shannon entropy'}
    )

fig.show()

#%% Correlation between scores

fig = px.scatter(score_df, x="Mutability", y="Shannon_entropy") #, trendline="ols")

fig.update_layout(
    title_text = 'Shannon entropy in function of Mutability per 39b for SARS-CoV-2',
    yaxis = {'title': 'Shannon entropy'},
    xaxis = {'title': 'Mutability'}
    )

fig.show()