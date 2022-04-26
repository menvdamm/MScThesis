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
if not os.path.isdir('./Data/Aligned/Human'): os.mkdir('./Data/Aligned/Human')
if not os.path.isdir('./Data/Tree'): os.mkdir('./Data/Tree')

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



# Time for Human_A_1:  min (FFT-large-NS-2 (Not tested.))
# Time for Human_A_2:  min (FFT-large-NS-2 (Not tested.))
# Time for Human_A_3:  min (FFT-large-NS-2 (Not tested.))?
# Time for Human_A_4:  min (FFT-large-NS-2 (Not tested.))?
# Time for Human_A_5:  min (FFT-large-NS-2 (Not tested.))
# Time for Human_A_6:  min (FFT-NS-2)
# Time for Human_A_7: 11.1 min (FFT-NS-2)
# Time for Human_A_8: 5.7 min (FFT-NS-2)

# Time for Human_B_1: 5.2 min (FFT-NS-2)
# Time for Human_B_2:  min (FFT-NS-2)
# Time for Human_B_3:  min (FFT-NS-2)
# Time for Human_B_4:  min (FFT-NS-2)
# Time for Human_B_5:  min (FFT-NS-2)
# Time for Human_B_6:  min (FFT-NS-2)
# Time for Human_B_7: 1.4 min (FFT-NS-2)
# Time for Human_B_8: 1.3 min (FFT-NS-2)


#%% tree output
# with open('Align_INFO.txt', 'w') as file:
#     for Type in ['A', 'B']:
#         for Segment in list('12345678'):
#             input_file = './Data/Clean/Human/'+Type+'_'+Segment+'.fasta'
#             output_file = './Data/Aligned/Human/'+Type+'_'+Segment+'.fasta'
#             cmd = 'mafft --thread -1 --auto --preservecase --treeout --parttree --reorder '+input_file+' > '+output_file
#             t2 = time.time()
#             subprocess.run(cmd, shell = True)
#             t3 = time.time()
#             print('Time for Human_'+Type+'_'+Segment+' alignment + parttree:', round((t3-t2)/60, 1), 'min')
#             minutes = str(round((t3-t2)/60, 1))
#             file.write('Time for Human_'+Type+'_'+Segment+' alignment + parttree: '+minutes+' min\n')


 #%% Fix FASTA headers

# # mafft truncates headers that are too long so we need to replace the headers
# # in the aligned files with the headers in the cleaned files

# for Type in ['A', 'B']:
#     for Segment in list('12345678'):
#         old_alignment = './Data_old/Aligned/Human/'+Type+'_'+Segment+'_temp.fasta'
#         clean = './Data/Clean/Human/'+Type+'_'+Segment+'.fasta'
#         new_alignment = './Data/Aligned/Human/'+Type+'_'+Segment+'.fasta'
#         with (open(old_alignment, 'r') as old_alignment, 
#               open(clean, 'r') as clean, 
#               open(new_alignment, 'w') as new_alignment):
#             al_recs = SeqIO.parse(old_alignment, 'fasta')
#             cl_recs = SeqIO.parse(clean, 'fasta')
#             for al_rec, cl_rec in zip(al_recs, cl_recs):
#                 new_rec = cl_rec
#                 new_rec.seq = al_rec.seq
#                 SeqIO.write(new_rec, new_alignment, 'fasta')
#        os.remove('./Data/Aligned/Human/'+Type+'_'+Segment+'_temp.fasta')
                