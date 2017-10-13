#!/usr/bin/python

# The script takes three arguments: file containing all PE reads, file containing all scaffolds, and file containing
# scaffolds that are >1kb.
# Created by Adi Lavy, 2017

import sys
import numpy
import csv
import subprocess
import argparse
from Bio import SeqIO
from collections import defaultdict


#get the arguments from the user
def main(**kwargs):
	pe_reads_file = kwargs["pereads"]
	assembly_file = kwargs["assmb"]
	min1000_file = kwargs["minimal"]
	
# Create dictionaries which will contain the results
	pe_reads = {}
	full_stats = {}
	min1000_stats = {}

###### start the calculations

	# go over the PE file and extract info
	lengths_PE = [len(seq_record.seq) for seq_record in SeqIO.parse(pe_reads_file, "fasta")]


	#total number of bp in raw reads
	pe_reads['total_bp'] = sum(lengths_PE)
	pe_reads['number'] = len(lengths_PE)
############################
	#go over the assembly file of full assembly and extract info
	lengths_full = [len(seq_record.seq) for seq_record in SeqIO.parse(assembly_file, "fasta")]
	lengths_full = sorted (lengths_full)

	#total bp in the assembly
	full_stats['total_bp'] = sum(lengths_full)
	#total number of scaffolds
	full_stats['number'] = len(lengths_full)
	#what is the longest contig
	full_stats['longest'] = lengths_full[-1]
	#average length of contigs
	full_stats['avg_len'] = numpy.mean(lengths_full)

############################
	#go over the assembly file of min1000 and extract info
	lengths_min1000 = [len(seq_record.seq) for seq_record in SeqIO.parse(min1000_file, "fasta")]
	lengths_min1000 = sorted (lengths_min1000)

	#total bp in the assembly
	min1000_stats['total_bp'] = sum(lengths_min1000)
	#total number of scaffolds
	min1000_stats['number'] = len(lengths_min1000)
	#what is the longest contig, no need as it is the same as in the full set of scaffolds
	#min1000_stats['longest'] = lengths_min1000[-1]
	#average length of contigs
	min1000_stats ['avg_len'] = numpy.mean(lengths_min1000)


###############################
#calcs for full assembly
	#percent of bp in full assembly / bp in PE
	full_stats['percent_bp_assembled'] = full_stats['total_bp'] / float(pe_reads['total_bp']* 100)
	
	#N50
	full_stats['N50'] = n50_calc(lengths_full)
	#calculate n50 from https://caramagnabosco.wordpress.com/2014/02/20/calculate-the-n50-of-assembled-contigs/


###############################
#calcs for min1000 file
	min1000_stats['percent_bp_assembled'] = min1000_stats['total_bp'] / float(pe_reads['total_bp'] * 100)


#scaffolds in min1000 / total scaffolds (%)
	min1000_stats['perc_scaffs_in1000_total'] = min1000_stats['number'] / float(full_stats['number'] * 100)

#N50
	min1000_stats['N50'] = n50_calc(lengths_min1000)

	headers=['Total bp','Number of reads','total bp in scaffolds','%bp assembled','# scaffolds','avg scaffold length','Longest scaffold','N50','>1kb scaff #bp','>1kb scaffolds %bp assembled','>1kb #scaffolds','>1kb scaff average length','>1kb scaff N50 ','%scaffolds >1kb']

	all_values=[]
	all_values.append(pe_reads['total_bp'])
	all_values.append(pe_reads['number'])
	all_values.append(full_stats['total_bp'])
	all_values.append(full_stats['percent_bp_assembled'])
	all_values.append(full_stats['number'])
	all_values.append(full_stats['avg_len'])
	all_values.append(full_stats['longest'])
	all_values.append(full_stats['N50'])
	all_values.append(min1000_stats['total_bp'])
	all_values.append(min1000_stats['percent_bp_assembled'])
	all_values.append(min1000_stats['number'])
	all_values.append(min1000_stats['avg_len'])
	all_values.append(min1000_stats['N50'])
	all_values.append(min1000_stats['perc_scaffs_in1000_total'])
	
	with open('assembly_stats.csv', 'w') as myfile:
		writer = csv.writer(myfile, dialect='excel')
		writer.writerow(headers)
		writer.writerow(all_values)
	myfile.close()

def n50_calc(scaff_file):
	#calculate n50 from https://caramagnabosco.wordpress.com/2014/02/20/calculate-the-n50-of-assembled-contigs/

	unique = []
	for entry in scaff_file:
		if not entry in unique:
			unique.append(entry)

	n50 = []

	for entry in unique:
		multiplier = scaff_file.count(entry) * entry
		for i in range(multiplier):
		  n50.append(entry)

	index = int(len(n50)/2)

	avg = []

	if index % 2==0:
		first = n50[index - 1]
		second = n50[index]
		avg.append(first)
		avg.append(second)
		n50 = numpy.mean (avg)
		return n50
	else:
		return n50[index -1]

#start the main body of the script
#get the arguments
if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="This script will return stats about the rawreads and assembled contigs. It takes three arguments: file containing all PE reads, file containing all scaffolds, and file containing scaffolds larger than 1kb")
	parser.add_argument("-p", "--pairedreads", type=str, required=True, help="File with PE fasta reads")
	parser.add_argument("-a", "--assembly", type=str, required=True, help="File with all scaffold in fasta format")
	parser.add_argument("-m", "--min1000", type=str, required=True, help="File with scaffolds >1kb in fasta format")
	args = parser.parse_args() 
	main(pereads=args.pairedreads, assmb=args.assembly, minimal=args.min1000)
	
