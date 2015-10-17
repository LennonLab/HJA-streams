###############################################################################
# Python Script to Trim Based on designated start and end
###############################################################################
# Written by Mario Muscarella
# Last Update 31 July 2013

# Directions:

from Bio import SeqIO
import sys
import glob

# change these numbers
start = 0 
end = 250

def trim_positions(records, start, end):
	for record in records:
		yield record[start:end]

files = glob.glob("*R1_001.fastq")

for x in files:
	original_seqs = SeqIO.parse(x, "fastq")
	trimmed_seqs = trim_positions(original_seqs, start, end)
	output_handle = open(x.replace(".fastq","")+".trim.fastq", "w")
	count = SeqIO.write(trimmed_seqs, output_handle, "fastq")
	output_handle.close()
	print "Trimmed %i reads in %s" % (count, x)

