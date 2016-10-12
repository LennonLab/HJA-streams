#!/bin/bash
#PBS -k o
#PBS -l nodes=1:ppn=8,vmem=50gb,walltime=02:00:00
#PBS -M wisnoski@indiana.edu
#PBS -m abe
#PBS -N mothur-summary.seqs
#PBS -j oe
cd /N/dc2/projects/Lennon_Sequences/2015_HJA-streams
module load gcc/4.9.2
module load mothur/1.36.1
mothur "#summary.seqs(fasta=hja_streams.trim.contigs.good.unique.good.filter.unique.precluster.fasta, accnos=hja_trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.accnos, processors=8)"
