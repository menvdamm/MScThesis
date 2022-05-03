#!/usr/bin/env bash

#PBS -l nodes=1:ppn=all
#PBS -l walltime=24:00:00
#PBS -l mem=100gb
#PBS -m abe

cd $HOME/data

module purge
module load ViennaRNA

RNAfold -p0 -d2 --noLP --noPS --infile=SARSCoV2_consensus.fasta --outfile=SARSCoV2.ifold