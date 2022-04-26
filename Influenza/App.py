#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Created on Sat Oct  2 12:17:07 2021

@author: mennovandamme
'''

#### Pipeline to create interactive app

#%% Dependencies

import json
import pandas as pd
# conda install -c plotly plotly
import plotly.express as px
# conda install dash
import dash
from dash.dependencies import Input, Output
import dash_core_components as dcc
import dash_html_components as html
#  installing dash-bio: https://stackoverflow.com/questions/41060382/using-pip-to-install-packages-to-anaconda-environment
# find anaconda path: https://docs.anaconda.com/anaconda/user-guide/tasks/integration/python-path/
import dash_bio as dashbio
    
#%% Loading the necessary data

for Type in ['A', 'B']:
    for Segment in list('12345678'): 
        with open('./Data/Dataframes/Constant/con_df_'+Type+'_'+Segment+'.csv') as file:
            globals()['con_df_'+Type+'_'+Segment] = pd.read_csv(file)
        with open('./Data/Dataframes/Variable/var_df_'+Type+'_'+Segment+'.csv') as file:
            globals()['var_df_'+Type+'_'+Segment] = pd.read_csv(file)  
        with open('./Data/Dataframes/Normalized/norm_var_df_'+Type+'_'+Segment+'.csv') as file:
            globals()['norm_var_df_'+Type+'_'+Segment] = pd.read_csv(file)
            
with open('./Data/Metadata/references.json') as file:
    references = json.load(file)

#%% Only Influenza: multiple df, multiple dropdowns, RNA structure

# things to fix: the animation takes too long, change frame time to make it like 5 secs in total
# axis should change together with the dropdown selections
# spaces between the dropdowns and the Graph
# nicer background color and overall lay out
# better explanation of the Graphs
# DON'T use forna container: it doesn't work well, slows the site down by a LOT, pic would be better prob

# data
for Type in ['A', 'B']:
    for Segment in list('12345678'): 
        with open('./Data/Dataframes/Constant/con_df_'+Type+'_'+Segment+'.csv') as file:
            globals()['con_df_'+Type+'_'+Segment] = pd.read_csv(file)
        with open('./Data/Dataframes/Variable/var_df_'+Type+'_'+Segment+'.csv') as file:
            globals()['var_df_'+Type+'_'+Segment] = pd.read_csv(file)  
        with open('./Data/Dataframes/Normalized/norm_var_df_'+Type+'_'+Segment+'.csv') as file:
            globals()['norm_var_df_'+Type+'_'+Segment] = pd.read_csv(file)

# dropdown options
Type_options = [{'label': Type, 'value': Type} for Type in list('AB')]
Segment_options = [{'label': Segment, 'value': Segment} for Segment in list('12345678')]
Data_options = [{'label': Data, 'value': Data} for Data in ['con', 'var', 'norm_var']]

# dash app setup
app = dash.Dash(__name__)

# layout
app.layout = html.Div([
    html.H1(
        children = ['Influenza'], style = {'text-align': 'center'}),
    html.Div([
        'Select a type',
        dcc.Dropdown(id = 'Type_dropdown', 
                     options = Type_options,
                     value = 'A',
                     clearable = False,
                     searchable = False,
                     optionHeight = 35,
                     style = {'width': '45%'}
                     )]),
    html.Div([
        'Select a genome segment',
        dcc.Dropdown(id = 'Segment_dropdown', 
                     options = Segment_options,
                     value = '1',
                     clearable = False,
                     searchable = False,
                     optionHeight = 35,
                     style = {'width': '45%'}
                     )]),
    html.Div([
        'Select which data to display',
        dcc.Dropdown(id = 'Data_dropdown', 
                     options = Data_options,
                     value = 'con',
                     clearable = False,
                     searchable = False,
                     optionHeight = 35,
                     style = {'width': '45%'}
                     )]),
    html.Div([dcc.Graph(id = 'Graph'),
    html.Br(),
    dashbio.FornaContainer(id = 'RNAfold',
                           height=800,
                           width=800)   
              ])
])

# callback
@app.callback(
    Output(component_id = 'Graph', component_property = 'figure'),
    Output(component_id = 'RNAfold', component_property = 'sequences'), 
    Input(component_id = 'Type_dropdown', component_property = 'value'),
    Input(component_id = 'Segment_dropdown', component_property = 'value'),
    Input(component_id = 'Data_dropdown', component_property = 'value')
    )

# Graph returning function
def display_Influenza_data(Type, Segment, Data):
        df = globals()[Data+'_df_'+Type+'_'+Segment]
        if Data == 'con':
            fig = px.bar(
                    data_frame = df,
                    x = 'Positions', y = 'Mutations', 
                    animation_frame = 'Year',
                    range_y = [0, max(df['Mutations']) + 100],
                    range_x = [1, max(df['Positions'])]
                    )
            fig.update_layout(
                    yaxis = {'title': 'Constant mutation count'},
                    xaxis = {'title': 'Position on the genome (bp)'}
                    )
    
        if Data == 'var':
            fig = px.bar(
                    data_frame = df,
                    x = 'Positions', y = 'Mutations', 
                    animation_frame = 'Year',
                    range_y = [0, max(df['Mutations']) + 100],
                    range_x = [1, max(df['Positions'])]
                    )
            fig.update_layout(
                    yaxis = {'title': 'Variable mutation count'},
                    xaxis = {'title': 'Position on the genome (bp)'}
                    ) 
        if Data == 'norm_var':
            fig = px.bar(
                    data_frame = df,
                    x = 'Positions', y = 'Mutations', 
                    animation_frame = 'Year',
                    range_y = [0, max(df['Mutations']) + 10],
                    range_x = [1, max(df['Positions'])]
                    )
            fig.update_layout(
                    yaxis = {'title': 'Normalized variable mutation count'},
                    xaxis = {'title': 'Position on the genome (bp)'}
                    )
        sequence = {'sequence': references[Type][Segment]['Sequence'],
                    'structure': references[Type][Segment]['Structure'],
                    'options': {'applyForce': False,
                                'avoidOthers': True,
                                'circulariseExternal': True,
                                'labelInterval': 50}}        
        return fig, [sequence]   

# run
if __name__ == '__main__':
    app.run_server(debug = False)    

#%% All viruses

# things to fix: the animation takes too long, change frame time to make it like 5 secs in total
# axis should change together with the dropdown selections
# spaces between the dropdowns and the Graph
# nicer background color and overall lay out
# better explanation of the Graphs
# change data type labels

# data
for Type in ['A', 'B']:
    for Segment in list('12345678'): 
        with open('./Data/Dataframes/Constant/con_df_'+Type+'_'+Segment+'.csv') as file:
            globals()['con_df_'+Type+'_'+Segment] = pd.read_csv(file)
        with open('./Data/Dataframes/Variable/var_df_'+Type+'_'+Segment+'.csv') as file:
            globals()['var_df_'+Type+'_'+Segment] = pd.read_csv(file)  
        with open('./Data/Dataframes/Normalized/norm_var_df_'+Type+'_'+Segment+'.csv') as file:
            globals()['norm_var_df_'+Type+'_'+Segment] = pd.read_csv(file)

with open('./Data/Metadata/references.json') as file:
    references = json.load(file)

# dropdown options
Type_options = [{'label': Type, 'value': Type} for Type in list('AB')]
Segment_options = [{'label': Segment, 'value': Segment} for Segment in list('12345678')]
Data_options = [{'label': 'Constant mutations', 'value': 'con'},
                {'label': 'Variable mutations', 'value': 'var'},
                {'label': 'Normalized variable mutations', 'value': 'norm_var'}]


# dash app setup
#external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css'] #doesn't work
app = dash.Dash(__name__)
#                external_stylesheets = external_stylesheets
                
# https://dash.plotly.com/callback-gotchas:
# app.config['suppress_callback_exceptions'] = True

# layout

Influenza_layout = html.Div(
    id = "Influenza_layout", 
    style = {'display': 'block'},
    children = [
        html.H3('Influenza'),
        html.Div([
            'Select a type',
            dcc.Dropdown(id = 'Type_dropdown',
                         options = Type_options,
                         value = 'A',
                         clearable = False,
                         searchable = False,
                         optionHeight = 35,
                         style = {'width': '45%'}
                         )]),
        html.Br(),
        html.Div([
            'Select a genome segment',
            dcc.Dropdown(id = 'Segment_dropdown', 
                         options = Segment_options,
                         value = '1',
                         clearable = False,
                         searchable = False,
                         optionHeight = 35,
                         style = {'width': '45%'}
                         )]),
        html.Br(),
        html.Div([
            'Select which data to display',
            dcc.Dropdown(id = 'Data_dropdown', 
                         options = Data_options,
                         value = 'con',
                         clearable = False,
                         searchable = False,
                         optionHeight = 35,
                         style = {'width': '45%'}
                         )]), 
            html.Br(),
            dcc.Graph(id = 'Graph'),
            html.Br(),
            dashbio.FornaContainer(id = 'RNAfold',
                                    height=800,
                                    width=800)                
            ])

Corona_layout = html.Div(
            id = 'Corona_layout',
            style = {'display': 'none'},
            children = [
                html.H3(children = ['Corona']),
                html.Div(['blablabla'])                
            ])

HIV_layout = html.Div(
    id = 'HIV_layout',
    style = {'display': 'none'},
    children = [
        html.H3(children = ['HIV']), 
        html.Div(['Select a type', 
        dcc.Dropdown(id = 'HIV_dropdown',
                                  options = Type_options,
                                  value = 'A',
                                  clearable = False,
                                  searchable = False,
                                  optionHeight = 35,
                                  style = {'width': '45%'})])
    ])        
        
    
app.layout = html.Div([
    html.H1(children = ['RNA virus evolution'], style = {'text-align': 'center'}),
    dcc.Tabs(id = 'Tabs',  
             children = [dcc.Tab(label = 'Influenza', value = 'Influenza'),
                         dcc.Tab(label = 'SARS-CoV-2', value = 'Corona'),
                         dcc.Tab(label = 'HIV 1', value = 'HIV'),
                         ],
             value = 'Influenza'
             ),
    # html.Div(id = 'Layout'),
    # dcc.Store(id = 'Tab', data = 'Influenza'),
    # dcc.Store(id = 'Type'),
    # dcc.Store(id = 'Segment'),
    # dcc.Store(id = 'Data')
    Influenza_layout,
    Corona_layout,
    HIV_layout
    ])


# callbacks

# @app.callback(
#     Output(component_id = 'Layout', component_property = 'children'),
#     Output(component_id = 'Type_dropdown', component_property = 'value'),
#     Output(component_id = 'Segment_dropdown', component_property = 'value'),
#     Output(component_id = 'Data_dropdown', component_property = 'value'),
#     Input(component_id = 'Tabs', component_property = 'value')
#     )
# def tabs_update_layout(Tab):
#     if Tab == 'Influenza':
#         return Influenza_layout, 'A', '1', 'con'
#     elif Tab == 'Corona':
#         return Corona_layout, 'A', '1', 'con'
#     elif Tab == 'HIV':
#         return HIV_layout, 'A', '1', 'con'

# @app.callback(
#     Output(component_id = 'Layout', component_property = 'children'),
#     Output(component_id = 'Tab', component_property = 'data'),
#     Input(component_id = 'Tabs', component_property = 'value')
#     )
# def tabs_update_layout(Tab):
#     if Tab == 'Influenza':
#         return Influenza_layout, Tab
#     elif Tab == 'Corona':
#         return Corona_layout, Tab
#     elif Tab == 'HIV':
#         return HIV_layout, Tab
#
#Influenza_layout['Type_dropdown'].value
#
# if app.layout['Tab'].data == 'Influenza':
    
    
@app.callback(
    Output(component_id = 'Influenza_layout', component_property = 'style'),
    Output(component_id = 'Corona_layout', component_property = 'style'),
    Output(component_id = 'HIV_layout', component_property = 'style'),
    Input(component_id = 'Tabs', component_property = 'value')
    )
def tabs_update_layout(Tab):
    if Tab == 'Influenza':
        return {'display': 'block'}, {'display': 'none'}, {'display': 'none'}
    elif Tab == 'Corona':
        return {'display': 'none'}, {'display': 'block'}, {'display': 'none'}
    elif Tab == 'HIV':
        return {'display': 'none'}, {'display': 'none'}, {'display': 'block'}
    
    
@app.callback(
        Output(component_id = 'Graph', component_property = 'figure'),
        Output(component_id = 'RNAfold', component_property = 'sequences'),    
        Input(component_id = 'Type_dropdown', component_property = 'value'),
        Input(component_id = 'Segment_dropdown', component_property = 'value'),
        Input(component_id = 'Data_dropdown', component_property = 'value')
        )
def display_Influenza_data(Type, Segment, Data):
        df = globals()[Data+'_df_'+Type+'_'+Segment]
        if Data == 'con':
            fig = px.bar(
                    data_frame = df,
                    x = 'Positions', y = 'Mutations', 
                    animation_frame = 'Year',
                    range_y = [0, max(df['Mutations']) + 100],
                    range_x = [1, max(df['Positions'])]
                    )
            fig.update_layout(
                    yaxis = {'title': 'Constant mutation count'},
                    xaxis = {'title': 'Position on the genome (bp)'}
                    )
    
        if Data == 'var':
            fig = px.bar(
                    data_frame = df,
                    x = 'Positions', y = 'Mutations', 
                    animation_frame = 'Year',
                    range_y = [0, max(df['Mutations']) + 100],
                    range_x = [1, max(df['Positions'])]
                    )
            fig.update_layout(
                    yaxis = {'title': 'Variable mutation count'},
                    xaxis = {'title': 'Position on the genome (bp)'}
                    ) 
        if Data == 'norm_var':
            fig = px.bar(
                    data_frame = df,
                    x = 'Positions', y = 'Mutations', 
                    animation_frame = 'Year',
                    range_y = [0, max(df['Mutations']) + 10],
                    range_x = [1, max(df['Positions'])]
                    )
            fig.update_layout(
                    yaxis = {'title': 'Normalized variable mutation count'},
                    xaxis = {'title': 'Position on the genome (bp)'}
                    )
        sequence = {'sequence': references[Type][Segment]['Sequence'],
                    'structure': references[Type][Segment]['Structure'],
                    'options': {'applyForce': False,
                                'avoidOthers': True,
                                'circulariseExternal': True,
                                'labelInterval': 50}}        
        return fig, [sequence]   
    
    
    
    # elif Tab == 'Corona':
    #     fig = 0
    #     sequence = 0
             
    # elif Tab == 'HIV':
    #     fig = 0
    #     sequence = 0
        
    # return fig, [sequence]

# run
if __name__ == '__main__':
    app.run_server(debug = False)
    

