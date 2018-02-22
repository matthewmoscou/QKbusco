# /usr/bin/python

# 1) Run BUSCO on a collection of genomes or transcriptomes
# 2) Assess the completeness of BUSCO analysis and individual gene families
#	Parameters: Complete/Duplicated/Fragmented; percent of the assemblies; how to handle duplicated (ignore, include)
# 3) Generate concatenated set of genes (locally aligned gene families - PRANK)
#	Transdecoder - identify ORFs
#	Select open reading frame
#	Merge gene families into individual FASTA

# open questions
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

## OptionParser
# import arguments and options
usage = "usage: %prog [BUSCO folders]"
parser = OptionParser(usage=usage)
parser.add_option("-b", "--busco", action="store", type="string", dest="busco", default='', help="BUSCO orthogroup information file")
parser.add_option("-m", "--master", action="store", type="string", dest="master", default='', help="Master file with identifier, BUSCO full table output, Transdecoder")
(options, args) = parser.parse_args()

## Functions
def parse_text_delimited(file_name, exclude, index):
	information = []
	Fopen = open(file_name, 'r')
	for line in Fopen.readlines():
		if len(line) > 0:
			if line[0] != exclude:
				sline = string.split(line)
				microinformation = []
				for indice in index:
					if indice < len(sline):
						microinformation.append(sline[indice])
				information.append(microinformation)
	Fopen.close()
	return information

def main():
	# inventory of individual BUSCO files
	BUSCO_inventory = parse_text_delimited(options.master, '#', range(5))
	BUSCO_gene_information = parse_text_delimited(options.busco, 'O', range(8))

	accession_BUSCO = {}

	# for every BUSCO analysis, identify genes
	for analysis in BUSCO_inventory:
		accession_BUSCO[analysis[0]] = [[], [], []]

		current_BUSCO = parse_text_delimited(analysis[2], '#', range(4))

		gene_duplicates = {}

		for gene in current_BUSCO:
			if len(gene) > 2:
				if gene[0] not in gene_duplicates.keys():
					gene_duplicates[gene[0]] = []

				gene_duplicates[gene[0]].append(gene)

		for gene in gene_duplicates.keys():
			if len(gene_duplicates[gene]) == 1:
				if gene_duplicates[gene][0][1] in ['Complete']:
					accession_BUSCO[analysis[0]][0].append(gene_duplicates[gene][0])
				if gene_duplicates[gene][0][1] in ['Fragmented']:
					accession_BUSCO[analysis[0]][1].append(gene_duplicates[gene][0])
			elif len(gene_duplicates[gene]) > 1:
				maxscore = 0
				maxgeneindex = -1

				for duplicate_index in range(len(gene_duplicates[gene])):
					if float(gene_duplicates[gene][duplicate_index][3]) > maxscore:
						maxgeneindex = duplicate_index

				if gene_duplicates[gene][maxgeneindex][1] == 'Duplicated':
					accession_BUSCO[analysis[0]][2].append(gene_duplicates[gene][maxgeneindex])
	
		print analysis[0], len(accession_BUSCO[analysis[0]][0]), len(accession_BUSCO[analysis[0]][1]), len(accession_BUSCO[analysis[0]][2])

	return
	
if __name__ == '__main__':
	main()
