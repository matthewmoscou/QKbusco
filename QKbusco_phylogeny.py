# /usr/bin/python

# 1) Generate concatenated set of genes (locally aligned gene families - PRANK)
#	PRANK (PHYLIP)
#	Concatenate aligned sequences
# 2) Phylogenetic analysis

# note
# the longest ORF may not always be the one of interest
#	if a major issue, we can parse the HMM output and identify the correct ORF, need to keep information from Transdecoder
# is there a case that two BUSCO genes hit the same transcriptome contig?
#	yes, this will need to be dealt with in the ORF prediction!
# are all fragmented genes single copy?

## Modules
import math

import optparse
from optparse import OptionParser 

import os.path, subprocess

import sets

import stats

import string

import time
from time import gmtime, strftime


## OptionParser
# import arguments and options
usage = "usage: %prog -d depth -s superalignment.phy [PHYLIP_output]"
parser = OptionParser(usage=usage)
parser.add_option("-d", "--depth", action="store", type="string", dest="depth", default='0.4', help="Depth required in alignments to include in master alignment")
parser.add_option("-s", "--superalignment", action="store", type="string", dest="superalignment", default='superalignment.phy', help="Super alignment for all species")
(options, args) = parser.parse_args()

## Functions
def convert_phylip(phylip_input, phylip_output):
	Fopen = open(phylip_input, 'r')
	Fout = open(phylip_output, 'w')

	line = Fopen.readline()

	genes, alignment_length = string.split(line)
	genes = int(genes)
	alignment_length = int(alignment_length)
	Fout.write(line)

	for gene_index in range(genes):
		line = Fopen.readline()
		gene = string.split(line)[0]
		species_alignment = ''

		for sequence_index in range(int(math.ceil(alignment_length / 60.0))):
			line = Fopen.readline()
			species_alignment += string.split(line)[0]

		Fout.write(gene)
		
		for index in range(10 - len(gene)):
			Fout.write(' ')

		Fout.write('  ')

		Fout.write(species_alignment + '\n')
	
	Fopen.close()
	Fout.close()
	return

def parse_phylip(phylip_input):
	species_alignment = {}

	Fopen = open(phylip_input, 'r')

	line = Fopen.readline()
	line = Fopen.readline()

	while line:
		species, alignment = string.split(line)

		species_alignment[species] = alignment

		line = Fopen.readline()
	
	Fopen.close()

	return species_alignment

def main():
	species_busco_genes = {}
	species_superalignment = {}

	# intialize species_busco_genes with species BUSCO hits
	for arg in args:
		print arg
		convert_phylip(arg, string.replace(arg, '.best.phy', ''))
		species_alignment = parse_phylip(string.replace(arg, '.best.phy', ''))

		for species in species_alignment.keys():
			species_information = string.split(species, '_')

			species_busco_genes[species_information[0]] = []

	# for every phylip file, convert to relaxed Phylip format
	# run QKphylogeny_alignment_analysis.py
	for arg in args:
		print 'superset:\t', arg, strftime("%Y-%m-%d %H:%M:%S", gmtime())

		process_name = "python QKphylogeny_alignment_analysis.py -q -a %s -d %s -o %s" % (string.replace(arg, '.best.phy', ''), options.depth, string.replace(arg, '.phy.best.phy', '.e.phy'))
		process = subprocess.Popen(process_name, shell = True)
		process.wait()

		# evaluate the super set of genes (for all species)
		species_alignment = parse_phylip(string.replace(arg, '.phy.best.phy', '.e.phy'))

		for species in species_alignment.keys():
			species_information = string.split(species, '_')

			species_busco_genes[species_information[0]].append(string.replace(arg, '.PRANK.phy.best.phy', ''))

	# concatenate all species, evaluate amount of missing data for every species
	print species_busco_genes.keys()
	for species in species_busco_genes.keys():
		species_superalignment[species] = ''

	for arg in args:
		print 'merge species:\t', arg, strftime("%Y-%m-%d %H:%M:%S", gmtime())
		# read multiple sequence alignment
		species_alignment = parse_phylip(string.replace(arg, '.phy.best.phy', '.e.phy'))

		# set alignment length (for missing data)
		alignment_length = len(species_alignment[species_alignment.keys()[0]])
		missing_template = ''

		for index in range(alignment_length):
			missing_template += '-'

		# add alignments for all species in alignment
		species_present = []

		for species in species_alignment.keys():
			species_information = string.split(species, '_')
			species_superalignment[species_information[0]] += species_alignment[species]
			species_present.append(species_information[0])

		# add missing data for all other species not in the alignment
		for species in list(sets.Set(species_superalignment.keys()) - sets.Set(species_present)):
			species_superalignment[species] += missing_template
	
	# print the degree of missing data for each species
	for species in species_superalignment.keys():
		print species, species_superalignment[species].count('-'), len(species_superalignment[species])

	# export concatenated alignment
	Fopen = open(options.superalignment, 'w')
	Fopen.write(' ' + str(len(species_superalignment.keys())) + ' ' + str(len(species_superalignment[species_superalignment.keys()[0]])) + '\n')
	
	for species in species_superalignment.keys():
		Fopen.write(species)

		for index in range(10 - len(species)):
			Fopen.write(' ')

		Fopen.write(species_superalignment[species] + '\n')

	Fopen.close()

	return
	
if __name__ == '__main__':
	main()

