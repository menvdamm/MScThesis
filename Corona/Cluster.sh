#!/usr/bin/env bash

#PBS -l nodes=1:ppn=8
#PBS -l walltime=24:00:00
#PBS -l mem=25gb
#PBS -m abe

cd $HOME/data

module purge
module load CD-HIT

cd-hit-est -i SARSCoV2.fasta -o SARSCoV2_9985.clstr -c 0.9985 -n 10 -T 0 -M 25000