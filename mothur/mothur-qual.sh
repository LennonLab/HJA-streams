#!/bin/bash
#PBS -k o
#PBS -l nodes=1:ppn=8,vmem=100gb,walltime=5:00:00
#PBS -M wisnoski@indiana.edu
#PBS -m abe
#PBS -N mothur-qualcheck
#PBS -j oe
cd /N/dc2/projects/Lennon_Sequences/2015_HJA-streams
module load gcc/4.9.2
module load mothur/1.36.1
mothur "#make.contigs(file=hja_streams.files, processors=8); summary.seqs(fasta=hja_streams.trim.contigs.fasta); screen.seqs(fasta=hja_streams.trim.contigs.fasta, group=hja_streams.contigs.groups, summary=hja_streams.trim.contigs.summary, maxambig=0, maxlength=275)"
