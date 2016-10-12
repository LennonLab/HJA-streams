#!/bin/bash
#PBS -k o
#PBS -l nodes=1:ppn=1,vmem=500gb,walltime=48:00:00
#PBS -M wisnoski@indiana.edu
#PBS -m abe
#PBS -j oe
cd /N/dc2/projects/Lennon_Sequences/2015_HJA-streams
module load gcc/4.9.2
module load boost/1.52.0
module load mothur/1.38.1
mothur make.tree.batch
