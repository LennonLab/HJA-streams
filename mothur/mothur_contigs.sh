#!/bin/bash
#PBS -k o
#PBS -l nodes=1:ppn=8,vmem=25gb,walltime=2:00:00
#PBS -M wisnoski@indiana.edu
#PBS -m abe
#PBS -N mothurcontigs
#PBS -j oe
cd /N/dc2/projects/Lennon_Sequences/2015_HJA-streams
module load gcc/4.9.2
module load mothur/1.36.1
mothur "#make.contigs(file=hja_streams.files, processors=8)"
