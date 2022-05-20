#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct  2 13:22:57 2021

@author: mennovandamme
"""

#### Creating alignments

#%% Dependencies

import time, os, subprocess
# conda install -c bioconda mafft

#%% Directories

if not os.path.isdir('./Data/Aligned'): os.mkdir('./Data/Aligned')

#%% Alignment

# the best mafft option for this amount of sequences per file and the length 
# of the sequences is probably this one:
#
# FFT-NS-1 (very fast but very rough; progressive method):
#   mafft --thread -1 --retree 1 --maxiterate 0 --preservecase input [> output]
#
# another option that is more accurate but takes about 2x as long
#
# FFT-NS-2 (fast; progressive method):
#   mafft --thread -1 --retree 2 --maxiterate 0 --preservecase input [> output]
#
# or the version that chooses a technique automatically:
#   mafft --thread -1 --auto --preservecase input_file > +output_file


# manual version
# Type = 'B'
# Segment = '7'
# input_file = './Data/Clean/Human/'+Type+'_'+Segment+'.fasta'
# output_file = './Data/Aligned/Human/'+Type+'_'+Segment+'_temp.fasta'
# cmd = 'mafft --thread -1 --auto --preservecase '+input_file+' > '+output_file
# t2 = time.time()
# subprocess.run(cmd, shell = True)
# t3 = time.time()
# print('Time for Human_'+Type+'_'+Segment+':', round((t3-t2)/60, 1), 'min')

with open('Align_INFO.txt', 'w') as file:
    for Type in ['A', 'B']:
        for Segment in list('12345678'):
            input_file = './Data/Clean/Human/'+Type+'_'+Segment+'.fasta'
            output_file = './Data/Aligned/Human/'+Type+'_'+Segment+'.fasta'
            cmd = 'mafft --thread -1 --auto --preservecase '+input_file+' > '+output_file
            t2 = time.time()
            subprocess.run(cmd, shell = True)
            t3 = time.time()
            print('Time for Human_'+Type+'_'+Segment+' alignment:', round((t3-t2)/60, 1), 'min')
            minutes = str(round((t3-t2)/60, 1))
            file.write('Time for Human_'+Type+'_'+Segment+' alignment: '+minutes+' min\n')