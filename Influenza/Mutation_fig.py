–#!/usr/bin/env python3
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
from plotly.subplots import make_subplots
# conda install statsmodels
import statsmodels
import plotly.io as pio
import numpy as np


# setting the renderer to browser
pio.renderers.default = 'browser'
pio.templates.default = "plotly_white"

#%% Variability per segment animation

for Type in ['A', 'B']:
    for Segment in list('12345678'): 
        with open('./Data/Dataframes/Animations/animation_df_'+Type+'_'+Segment+'.csv') as file:
            globals()['animation_df_'+Type+'_'+Segment] = pd.read_csv(file)

Type = 'B'
Segment = '8'

Color = ['#000000']*len(globals()['animation_df_'+Type+'_'+Segment]['Position'])

# making animated Graph
fig = px.bar(data_frame = globals()['animation_df_'+Type+'_'+Segment], 
             x = 'Position', y = 'Variability', 
             animation_frame = 'Year',
             animation_group = 'Position',
             range_y = [0, max(globals()['animation_df_'+Type+'_'+Segment]['Variability']) + 2],
             color = Color,
             color_discrete_map = "identity"
             )

fig.update_traces(marker_color = 'black', 
                  marker_line_color = 'black',
                  marker_line_width = 0.05)

fig.update_layout(transition = {'duration': 10})

# changing axis titles
fig.update_layout(
    title_text = 'Influenza '+Type+' genome segment '+Segment,
    yaxis = {'title': 'Variability'},
    xaxis = {'title': 'Position'}
    )

# showing the Graph in browser
fig.show()

#%% Variability per type animation

## PROBLEM: animation is very slow, might be fixable by using GO figures from plotly

for Type in ['A', 'B']:
    with open('./Data/Dataframes/Animations/animation_df_'+Type+'.csv') as file:
        globals()['animation_df_'+Type] = pd.read_csv(file)

Type = 'A'

Color = ['#000000']*len(globals()['animation_df_'+Type]['Position'])

# making animated Graph
fig = px.bar(data_frame = globals()['animation_df_'+Type], 
             x = 'Position', y = 'Variability', 
             animation_frame = 'Year',
             animation_group = 'Position',
             range_y = [0, max(globals()['animation_df_'+Type]['Variability']) + 2],
             color = Color,
             color_discrete_map = "identity"
             )

fig.update_traces(marker_color = 'black', 
                  marker_line_color = 'black',
                  marker_line_width = 0.05)

fig.update_layout(transition = {'duration': 50})

# changing axis titles
fig.update_layout(
    title_text = 'Influenza '+Type,
    yaxis = {'title': 'Variability'},
    xaxis = {'title': 'Position'}
    )

# showing the Graph in browser
fig.show()

#%% Variability per type bar chart

for Type in ['A', 'B']:
    with open('./Data/Dataframes/Ungapped/df_'+Type+'.csv') as file:
        globals()['df_'+Type] = pd.read_csv(file) 
        
Type = 'A'

Color = ['#000000']*len(globals()['df_'+Type]['Position'])

# making animated Graph
fig = px.bar(data_frame = globals()['df_'+Type], 
             x = 'Position', y = 'Variability',
             range_y = [0, max(globals()['df_'+Type]['Variability']) + 2],
             color = Color,
             color_discrete_map = "identity"
             )

fig.update_traces(marker_color = 'black', 
                  marker_line_color = 'black',
                  marker_line_width = 0.05)

# changing axis titles
fig.update_layout(
    title_text = 'Variability per position for Influenza '+Type,
    yaxis = {'title': 'Variability'},
    xaxis = {'title': 'Position'}
    )

# showing the Graph in browser
fig.show()

#%% Mutation count per type bar chart

for Type in ['A', 'B']:
    with open('./Data/Dataframes/Ungapped/df_'+Type+'.csv') as file:
        globals()['df_'+Type] = pd.read_csv(file) 
        
Type = 'A'

Color = ['#000000']*len(globals()['df_'+Type]['Position'])

# making animated Graph
fig = px.bar(data_frame = globals()['df_'+Type], 
             x = 'Position', y = 'Mutation_count',
             range_y = [0, max(globals()['df_'+Type]['Mutation_count']) + 2],
             color = Color,
             color_discrete_map = "identity"
             )

fig.update_traces(marker_color = 'black', 
                  marker_line_color = 'black',
                  marker_line_width = 0.05)

# changing axis titles
fig.update_layout(
    title_text = 'Mutation count per position for Influenza '+Type,
    yaxis = {'title': 'Mutation count'},
    xaxis = {'title': 'Position'}
    )

# showing the Graph in browser
fig.show()

#%% Shannon entropy per type bar chart

for Type in ['A', 'B']:
    with open('./Data/Dataframes/Ungapped/df_'+Type+'.csv') as file:
        globals()['df_'+Type] = pd.read_csv(file) 
        
Type = 'A'

Color = ['#000000']*len(globals()['df_'+Type]['Position'])

# making animated Graph
fig = px.bar(data_frame = globals()['df_'+Type], 
             x = 'Position', y = 'Shannon_entropy',
             range_y = [0, max(globals()['df_'+Type]['Shannon_entropy']) + 0.5],
             color = Color,
             color_discrete_map = "identity"
             )

fig.update_traces(marker_color = 'black', 
                  marker_line_color = 'black',
                  marker_line_width = 0.05)

# changing axis titles
fig.update_layout(
    title_text = 'Shannon entropy per position for Influenza '+Type,
    yaxis = {'title': 'Shannon entropy'},
    xaxis = {'title': 'Position'}
    )

# showing the Graph in browser
fig.show()

#%% Variability histogram

for Type in ['A', 'B']:
    with open('./Data/Dataframes/Ungapped/df_'+Type+'.csv') as file:
        globals()['df_'+Type] = pd.read_csv(file) 

Type = 'A'

fig = px.histogram(globals()['df_'+Type]['Variability'], nbins = 300)

fig.update_layout(
    title_text = 'Variability histogram for Influenza '+Type,
    yaxis = {'title': 'Position count'},
    xaxis = {'title': 'Variability'}
    )

fig.show()

#%% Mutation count histogram

for Type in ['A', 'B']:
    with open('./Data/Dataframes/Ungapped/df_'+Type+'.csv') as file:
        globals()['df_'+Type] = pd.read_csv(file) 

Type = 'A'
    
fig = px.histogram(globals()['df_'+Type]['Mutation_count'], nbins = 300)

fig.update_layout(
    title_text = 'Mutation count histogram for Influenza '+Type,
    yaxis = {'title': 'Position count'},
    xaxis = {'title': 'Mutation count'}
    )


fig.show()

#%% Shannon entropy histogram

for Type in ['A', 'B']:
    with open('./Data/Dataframes/Ungapped/df_'+Type+'.csv') as file:
        globals()['df_'+Type] = pd.read_csv(file) 

Type = 'A'
    
fig = px.histogram(globals()['df_'+Type]['Shannon_entropy'], nbins = 300)

fig.update_layout(
    title_text = 'Shannon entropy histogram for Influenza '+Type,
    yaxis = {'title': 'Position count'},
    xaxis = {'title': 'Shannon entropy'}
    )


fig.show()

#%% Correlation between scores

for Type in ['A', 'B']:
    with open('./Data/Dataframes/Ungapped/df_'+Type+'.csv') as file:
        globals()['df_'+Type] = pd.read_csv(file) 

Type = 'A'

fig = px.scatter(globals()['df_'+Type], x="Mutation_count", y="Variability") #, trendline="ols")

fig.update_layout(
    title_text = 'Mutation count in function of Variability per position for Influenza '+Type,
    yaxis = {'title': 'Shannon entropy'},
    xaxis = {'title': 'Mutation count'}
    )

fig.show()

fig = px.scatter(globals()['df_'+Type], x="Mutation_count", y="Shannon_entropy") #, trendline="ols")

fig.update_layout(
    title_text = 'Shannon entropy in function of Mutation count per position for Influenza '+Type,
    yaxis = {'title': 'Shannon entropy'},
    xaxis = {'title': 'Mutation count'}
    )

fig.show()

fig = px.scatter(globals()['df_'+Type], x="Variability", y="Shannon_entropy") #, trendline="ols")

fig.update_layout(
    title_text = 'Shannon entropy in function of Variability per position for Influenza '+Type,
    yaxis = {'title': 'Shannon entropy'},
    xaxis = {'title': 'Variability'}
    )

fig.show()

# fig.show()

#%% Shannon entropy bar chart extended

for Type in ['A', 'B']:
    with open('./Data/Dataframes/Ungapped/df_'+Type+'.csv') as file:
        globals()['df_'+Type] = pd.read_csv(file) 
    with open('./Data/Dataframes/Scores/scores_df_'+Type+'.csv') as file:
        globals()['scores_df_'+Type+'_39'] = pd.read_csv(file)
    for Segment in list('12345678'): 
         with open('./Data/Dataframes/Ungapped/df_'+Type+'_'+Segment+'.csv') as file:
             globals()['df_'+Type+'_'+Segment] = pd.read_csv(file)   

Type = 'B'

fig = go.Figure()

fig.add_trace(
    go.Bar(
        x = globals()['df_'+Type]['Position'],
        y = globals()['df_'+Type]['Shannon_entropy'],
        marker_color = 'black',
        marker_line_color = 'black',
        marker_line_width = 0.05,
        showlegend = False))

# for var in ['Variability_score', 'Mutation_count_score', 'Shannon_entropy_score']:
#     df = globals()['scores_df_'+Type+'_21'].sort_values(var)[0:10]
#     for index, row in df.iterrows():
#         fig.add_trace(
#             go.Scatter(
#                 x = [row['Begin_position'], row['End_position']],
#                 y = [-0.4, -0.4],
#                 mode = 'lines',
# #                name = '21 bp',
#                 line = dict(
#                     color = 'red', 
#                     width = 6)))
        
for var in ['Variability_score', 'Mutation_count_score', 'Shannon_entropy_score']:
    df = globals()['scores_df_'+Type+'_39'].sort_values(var)[0:10]
    for index, row in df.iterrows():
        fig.add_trace(
            go.Scatter(
                x = [row['Begin_position'], row['End_position']],
                y = [-0.1, -0.1],
                mode = 'lines',
#                name = '39 bp',
                line = dict(
                    color = 'blue', 
                    width = 6)))
        
# references['A']['1']['Features']

fig.add_trace(
    go.Scatter(
        x = [27, 2307],
        y = [-0.2, -0.2],
        mode = 'lines',
        line = dict(
            color = 'green', 
            width = 10),
        text = 'PB2',
        textposition = "middle center",
        textfont = dict(
            family="sans serif",
            size=18,
            color="black")))

fig.add_trace(
    go.Scatter(
        x = [24+len(df_A_1), 2298+len(df_A_1)],
        y = [-0.2, -0.2],
        mode = 'lines',
        line = dict(
            color = 'green', 
            width = 10),
        text = 'PB1',
        textposition = "middle center",
        textfont = dict(
            family="sans serif",
            size=18,
            color="black")))

fig.add_trace(
    go.Scatter(
        x = [24+len(df_A_1)+len(df_A_2), 2175+len(df_A_1)+len(df_A_2)],
        y = [-0.2, -0.2],
        mode = 'lines',
        line = dict(
            color = 'green', 
            width = 10),
        text = 'PA',
        textposition = "middle center",
        textfont = dict(
            family="sans serif",
            size=18,
            color="black")))

fig.add_trace(
    go.Scatter(
        x = [29+len(df_A_1)+len(df_A_2)+len(df_A_3), 1730+len(df_A_1)+len(df_A_2)+len(df_A_3)],
        y = [-0.2, -0.2],
        mode = 'lines',
        line = dict(
            color = 'green', 
            width = 10),
        text = 'HA',
        textposition = "middle center",
        textfont = dict(
            family="sans serif",
            size=18,
            color="black")))

length4 = len(df_A_1)+len(df_A_2)+len(df_A_3)+len(df_A_4)

fig.add_trace(
    go.Scatter(
        x = [45+length4, 1542+length4],
        y = [-0.2, -0.2],
        mode = 'lines',
        line = dict(
            color = 'green', 
            width = 10),
        text = 'NP',
        textposition = "middle center",
        textfont = dict(
            family="sans serif",
            size=18,
            color="black")))

length5 = len(df_A_1)+len(df_A_2)+len(df_A_3)+len(df_A_4)+len(df_A_5)

fig.add_trace(
    go.Scatter(
        x = [19+length5, 1429+length5],
        y = [-0.2, -0.2],
        mode = 'lines',
        line = dict(
            color = 'green', 
            width = 10),
        text = 'NA',
        textposition = "middle center",
        textfont = dict(
            family="sans serif",
            size=18,
            color="black")))

length6 = len(df_A_1)+len(df_A_2)+len(df_A_3)+len(df_A_4)+len(df_A_5)+len(df_A_6)

fig.add_trace(
    go.Scatter(
        x = [25+length6, 51+length6],
        y = [-0.2, -0.2],
        mode = 'lines',
        line = dict(
            color = 'green', 
            width = 10),
        text = 'M2',
        textposition = "middle center",
        textfont = dict(
            family="sans serif",
            size=18,
            color="black")))

fig.add_trace(
    go.Scatter(
        x = [739+length6, 1007+length6],
        y = [-0.2, -0.2],
        mode = 'lines',
        line = dict(
            color = 'green', 
            width = 10),
        text = 'M2',
        textposition = "middle center",
        textfont = dict(
            family="sans serif",
            size=18,
            color="black")))

fig.add_trace(
    go.Scatter(
        x = [25+length6, 784+length6],
        y = [-0.2, -0.2],
        mode = 'lines',
        line = dict(
            color = 'green', 
            width = 10),
        text = 'M1',
        textposition = "middle center",
        textfont = dict(
            family="sans serif",
            size=18,
            color="black")))

length7 = len(df_A_1)+len(df_A_2)+len(df_A_3)+len(df_A_4)+len(df_A_5)+len(df_A_6)+len(df_A_7)

fig.add_trace(
    go.Scatter(
        x = [25+length7, 56+length7],
        y = [-0.2, -0.2],
        mode = 'lines',
        line = dict(
            color = 'green', 
            width = 10),
        text = 'NS2',
        textposition = "middle center",
        textfont = dict(
            family="sans serif",
            size=18,
            color="black")))

fig.add_trace(
    go.Scatter(
        x = [528+length7, 864+length7],
        y = [-0.2, -0.2],
        mode = 'lines',
        line = dict(
            color = 'green', 
            width = 10),
        text = 'NS2',
        textposition = "middle center",
        textfont = dict(
            family="sans serif",
            size=18,
            color="black")))

fig.add_trace(
    go.Scatter(
        x = [26+length7, 719+length7],
        y = [-0.2, -0.2],
        mode = 'lines',
        line = dict(
            color = 'green', 
            width = 10),
        text = 'NS1',
        textposition = "middle center",
        textfont = dict(
            family="sans serif",
            size=18,
            color="black")))

fig.update_layout(title_text = 'Shannon entropy per position for Influenza '+Type,
                  showlegend = False)

fig.update_yaxes(title_text = "Shannon entropy")
fig.update_xaxes(title_text = "Position")

fig.show()

