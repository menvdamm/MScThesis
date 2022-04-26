#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  2 13:43:46 2021

@author: mennovandamme
"""

### Finding the lowest variability regions

#%% Dependencies

import json
import pandas as pd
import numpy as np

#%% Loading the necessary data

with open('./Data/Metadata/references.json') as file:
    references = json.load(file)
    
for Type in ['A', 'B']:
    with open('./Data/Dataframes/Ungapped/df_'+Type+'.csv') as file:
        globals()['df_'+Type] = pd.read_csv(file)
    for Segment in list('12345678'):
        with open('./Data/Dataframes/Ungapped/df_'+Type+'_'+Segment+'.csv') as file:
            globals()['df_'+Type+'_'+Segment] = pd.read_csv(file)
        
#%% Score dataframes

def scores_df_maker(Type, size):
    score_df = pd.DataFrame({'Segment': [],
                             'Segment_begin_position': [],
                             'Segment_end_position': [],
                             'Begin_position': [],
                             'End_position': [],
                             'Variability_score': [],
                             'Mutation_count_score': [],
                             'Shannon_entropy_score': [],
                             'Gap_percentage': []})
    length = 0
    for Segment in list('12345678'):
        df = globals()['df_'+Type+'_'+Segment]
        seq = references[Type][Segment]['Ungapped_sequence']
        Begin = pd.Series(list(range(0, len(seq)-size+1)))
        End = pd.Series(list(range(size-1, len(seq))))
        Seg = [int(Segment)] * len(Begin)
        Var_scores = pd.Series([sum(df['Variability'][begin:end+1]) for begin, end in zip(Begin, End)])
        Mut_scores = pd.Series([sum(df['Mutation_count'][begin:end+1]) for begin, end in zip(Begin, End)])
        Ent_scores = pd.Series([sum(df['Shannon_entropy'][begin:end+1]) for begin, end in zip(Begin, End)])
        Gap_scores = pd.Series([sum(df['Gap_percentage'][begin:end+1]) for begin, end in zip(Begin, End)])
        segment_df = pd.DataFrame({'Segment': Seg,
                                   'Segment_begin_position': Begin + 1,
                                   'Segment_end_position': End + 1,
                                   'Begin_position': Begin + 1 + length,
                                   'End_position': End + 1 + length,     
                                   'Variability_score': Var_scores,
                                   'Mutation_count_score': Mut_scores,
                                   'Shannon_entropy_score': Ent_scores,
                                   'Gap_percentage': Gap_scores})
        length += len(seq)
        score_df = score_df.append(segment_df)
    score_df = score_df.sort_values('Variability_score')
    return score_df

for Type in ['A', 'B']:
    globals()['scores_df_'+Type+'_39'] = scores_df_maker(Type, size = 39)
    globals()['scores_df_'+Type+'_39'].to_csv('./Data/Dataframes/Scores/scores_df_'+Type+'.csv', index = False)
 
#%% Histograms of scores

for Type in ['A', 'B']:
    with open('./Data/Dataframes/Scores/scores_df_'+Type+'.csv') as file:
        globals()['scores_df_'+Type] = pd.read_csv(file) 
        
Type = 'A'
size = '39'

import plotly.express as px
import plotly.io as pio

pio.renderers.default = 'browser'
pio.templates.default = "plotly_white"


fig = px.histogram(globals()['scores_df_'+Type]['Variability_score'], nbins = 300)
fig.update_layout(
    title_text = 'Variability histogram for Influenza '+Type+' per '+size+' bp',
    yaxis = {'title': 'Count'},
    xaxis = {'title': 'Variability score'}
    )
fig.show()

fig = px.histogram(globals()['scores_df_'+Type]['Mutation_count_score'], nbins = 300)
fig.update_layout(
    title_text = 'Mutation count histogram for Influenza '+Type+' per '+size+' bp',
    yaxis = {'title': 'Count'},
    xaxis = {'title': 'Mutation count score'}
    )
fig.show()

fig = px.histogram(globals()['scores_df_'+Type]['Shannon_entropy_score'], nbins = 300)
fig.update_layout(
    title_text = 'Shannon entropy histogram for Influenza '+Type+' per '+size+' bp',
    yaxis = {'title': 'Count'},
    xaxis = {'title': 'Shannon entropy score'}
    )
fig.show()

import statsmodels

fig = px.scatter(globals()['scores_df_'+Type], x="Mutation_count_score", y="Variability_score", trendline="ols")

fig.show()

fig = px.scatter(globals()['scores_df_'+Type], x="Mutation_count_score", y="Shannon_entropy_score", trendline="ols")

fig.show()

fig = px.scatter(globals()['scores_df_'+Type], x="Shannon_entropy_score", y="Variability_score", trendline="ols")

fig.show()

#%% Finding the lowest scoring sequences

scores_df_A.head(10)
minimum = scores_df_A.loc[scores_df_A['Variability_score'] == min(scores_df_A['Variability_score'])]
begin = int(minimum['Begin_position'].iloc[0])
end = int(minimum['End_position'].iloc[0])
RNA = references['A']['1']['Ungapped_sequence'][begin-1:end]
no_gap_df_A.loc[begin-1:end-1]


scores_df_B.head(20)
minimum = scores_df_B.loc[scores_df_B['Variability_score'] == min(scores_df_B['Variability_score'])]
begin = minimum['Begin_position'].iloc[0]
end = minimum['End_position'].iloc[0]
RNA = ref_B[begin-1:end]
no_gap_df_B.loc[begin-1:end-1]
RNA
'TTTGGAGACACAATTGCCTAC' in references['B']['7']['Sequence']

