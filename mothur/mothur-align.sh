#!/bin/bash
#PBS -k o
#PBS -l nodes=1:ppn=8,vmem=250gb,walltime=48:00:00
#PBS -M wisnoski@indiana.edu
#PBS -m abe
#PBS -N mothur-align.seqs
#PBS -j oe
cd /N/dc2/projects/Lennon_Sequences/2015_HJA-streams
module load gcc/4.9.2
module load mothur/1.36.1
mothur "#align.seqs(fasta=hja_streams.trim.contigs.good.unique.fasta, reference=silva.v4.fasta, processors=8)"
