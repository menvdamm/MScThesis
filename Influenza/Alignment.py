#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct  2 13:22:57 2021

@author: mennovandamme
"""

#### Creating alignments

#%% Dependencies

import os, subprocess
# conda install -c bioconda mafft

#%% Directories

if not os.path.isdir('./Data/Aligned'): os.mkdir('./Data/Aligned')

#%% Alignment

for Type in ['A', 'B']:
    for Segment in list('12345678'):
        input_file = './Data/Clean/'+Type+'_'+Segment+'.fasta'
        output_file = './Data/Aligned/'+Type+'_'+Segment+'.fasta'
        cmd = 'mafft --thread -1 --auto --preservecase '+input_file+' > '+output_file
        subprocess.run(cmd, shell = True)