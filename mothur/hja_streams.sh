#!/bin/bash
#PBS -k o
#PBS -l nodes=1:ppn=8,vmem=250gb,walltime=168:00:00
#PBS -M wisnoski@indiana.edu
#PBS -m abe
#PBS -j oe
cd /N/dc2/projects/Lennon_Sequences/hja-streams
module load gcc/6.3.0
/N/u/wisnoski/Carbonate/mothur/mothur /gpfs/home/w/i/wisnoski/Carbonate/GitHub/HJA-streams/mothur/hja_streams.batch
