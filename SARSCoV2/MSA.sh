#!/usr/bin/env bash

#PBS -l nodes=3:ppn=all
#PBS -l walltime=72:00:00
#PBS -l mem=150gb
#PBS -m abe

cd data

module purge
module load MAFFT

mafft --thread -1 --auto --preservecase --6merpair --addfragments SARSCoV2_rest.fasta SARSCoV2_ref.fasta > SARSCoV2.fasta