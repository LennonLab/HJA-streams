#!/bin/bash 
#PBS -k o
#PBS -l nodes=1:ppn=1,vmem=32gb
#PBS -l walltime=120:00:00
#PBS -m abe
#PBS -M wisnoski@indiana.edu
#PBS -N hjanullmodels
#PBS -j oe
#PBS -V

cd /N/u/wisnoski/Karst/GitHub/HJA-streams

module load nlopt
module load curl
module load java/1.8.0_74
module load r

R CMD BATCH analysis/PhylogeneticEcology.R
