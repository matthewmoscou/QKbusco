# /usr/bin/python

# 1) Run BUSCO on a collection of genomes or transcriptomes
# 2) Assess the completeness of BUSCO analysis and individual gene families
#	Parameters: Complete/Duplicated/Fragmented; percent of the assemblies; how to handle duplicated (ignore, include)
# 3) Generate concatenated set of genes (locally aligned gene families - PRANK)
#	Transdecoder - identify ORFs
#	Select open reading frame
#	Merge gene families into individual FASTA
#	PRANK (PHYLIP)
#	Concatenate aligned sequences
# 4) Phylogenetic analysis

# note
# the longest ORF may not always be the one of interest
#	if a major issue, we can parse the HMM output and identify the correct ORF, need to keep information from Transdecoder
# is there a case that two BUSCO genes hit the same transcriptome contig?
#	yes, this will need to be dealt with in the ORF prediction!
# are all fragmented genes single copy?

## Modules
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
parser.add_option("-s", "--status", action="store", type="string", dest="status", default='', help="Strict or relaxed assessment using complete, duplicated, or fragmented genes, options for parameter include: complete or allcomplete or allfragmented")
parser.add_option("-o", "--output", action="store", type="string", dest="output", default='1', help="Output folder")
parser.add_option("-t", "--threshold", action="store", type="int", dest="threshold", default=2, help="Number of species required to evaluate a BUSCO gene")
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

def parse_Transdecoder(FASTA_file_name):
	FASTA_file = open(FASTA_file_name, 'r')
	FASTA_data = FASTA_file.readlines()
	gene_model_sequence = {}

	# initialize dictionary
	for line in FASTA_data:
		if len(line) > 0:
			if line[0] == '>':
				identifiers = string.split(line[1:], ':')
				ID = identifiers[2]
				gene_model_sequence[ID] = {}
	FASTA_file.close()

	# import sequence
	for line in FASTA_data:
		sline = string.split(line)
		if len(line) > 0:
			if line[0] == '>':
				identifiers = string.split(line[1:], ':')
				ID = identifiers[2]
				model = string.replace(identifiers[6], ' type', '')
				position = identifiers[9]
				start_position = int(position[:position.index('-')])
				stop_position = int(position[(position.index('-') + 1):position.index('(')])
				strand = position[position.index('(') + 1]

				# review later to confirm it is functioning
				if strand == '+':
					frame = min([start_position, stop_position]) % 3 - 1
				elif strand == '-':
					frame = min([start_position, stop_position]) % 3 + 4
					
				gene_model_sequence[ID][model] = [frame, '']
			else:
				gene_model_sequence[ID][model][1] += sline[0]
	FASTA_file.close()
	return gene_model_sequence

def main():
	# inventory of individual BUSCO files
	BUSCO_inventory = parse_text_delimited(options.master, '#', range(4))

	# all BUSCO genes and biological information 
	BUSCO_gene_information = parse_text_delimited(options.busco, 'O', range(8))
	BUSCO_gene_species = {}

	for information in BUSCO_gene_information:
		BUSCO_gene_species[information[0]] = []

	# for every BUSCO analysis, identify genes
	accession_BUSCO = {}
	accession_BUSCO_transcript = {}
	duplicated_BUSCO = []

	for analysis in BUSCO_inventory:
		accession_BUSCO[analysis[0]] = []
		accession_BUSCO_transcript[analysis[0]] = {}

		current_BUSCO = parse_text_delimited(analysis[1], '#', range(4))

		if options.status in ['complete', 'allfragmented']:
			for gene in current_BUSCO:
				if gene[1] == 'Complete':
					accession_BUSCO[analysis[0]].append(gene)
				if options.status == 'allfragmented':
					if gene[1] == 'Fragmented':
						accession_BUSCO[analysis[0]].append(gene)

		if options.status in ['allcomplete', 'allfragmented']:
			gene_duplicates = {}

			for gene in current_BUSCO:
				if len(gene) > 2:
					if gene[0] not in gene_duplicates.keys():
						gene_duplicates[gene[0]] = []

					gene_duplicates[gene[0]].append(gene)

			for gene in gene_duplicates.keys():
				if len(gene_duplicates[gene]) == 1:
					if options.status == 'allcomplete':
						if gene_duplicates[gene][0][1] in ['Complete']:
							accession_BUSCO[analysis[0]].append(gene_duplicates[gene][0])
					else:
						if gene_duplicates[gene][0][1] in ['Complete', 'Fragmented']:
							accession_BUSCO[analysis[0]].append(gene_duplicates[gene][0])
				elif len(gene_duplicates[gene]) > 1:
					maxscore = 0
					maxgeneindex = -1

					for duplicate_index in range(len(gene_duplicates[gene])):
						if float(gene_duplicates[gene][duplicate_index][3]) > maxscore:
							maxgeneindex = duplicate_index

					if gene_duplicates[gene][maxgeneindex][1] == 'Duplicated':
						accession_BUSCO[analysis[0]].append(gene_duplicates[gene][maxgeneindex])
		
		if options.status not in ['complete', 'allcomplete', 'allfragmented']:
			print 'Error - Incorrect status'
			return
	
		# identify trancripts matching multiple BUSCO genes
		transcript_BUSCO_gene = {}

		for gene in current_BUSCO:
			if len(gene) > 2:
				if gene[2] not in transcript_BUSCO_gene.keys():
					transcript_BUSCO_gene[gene[2]] = []
				transcript_BUSCO_gene[gene[2]].append(gene[0])

		for transcript in transcript_BUSCO_gene.keys():
			if len(transcript_BUSCO_gene[transcript]) > 1:
				#print transcript, transcript_BUSCO_gene[transcript]
				for BUSCO_gene in transcript_BUSCO_gene[transcript]:
					duplicated_BUSCO.append(BUSCO_gene)

		for gene in accession_BUSCO[analysis[0]]:
			BUSCO_gene_species[gene[0]].append(analysis[0])
			accession_BUSCO_transcript[analysis[0]][gene[0]] = gene[2]
	
	gene_species_representation = []
	BUSCO_genes_for_analysis = []

	for gene in BUSCO_gene_species.keys():
		gene_species_representation.append(len(sets.Set(BUSCO_gene_species[gene])))

		if len(sets.Set(BUSCO_gene_species[gene])) >= options.threshold:
			BUSCO_genes_for_analysis.append(gene)

	BUSCO_genes_for_analysis = list(sets.Set(BUSCO_genes_for_analysis) - sets.Set(duplicated_BUSCO))

	print 'Minimum number of species observed for BUSCO genes:', min(gene_species_representation)
	print 'Maximum number of species observed for BUSCO genes:', max(gene_species_representation)
	print 'Median number of species observed for BUSCO genes:', stats.median(gene_species_representation)
	print 'Mean number of species observed for BUSCO genes:', stats.mean(gene_species_representation)
	print 'Number of genes with maximum number of species versus total number of BUSCO genes:', gene_species_representation.count(max(gene_species_representation)), '/', len(gene_species_representation)
	print 'Number of BUSCO genes with ambiguous mapping:', len(sets.Set(duplicated_BUSCO))

	# import open reading frames
	accession_gene_ORF = {}
	short_BUSCO = []

	for analysis in BUSCO_inventory:
		accession_gene_ORF[analysis[0]] = {}
		all_ORFs = parse_Transdecoder(analysis[2])

		for gene in accession_BUSCO[analysis[0]]:
			if gene[0] not in duplicated_BUSCO:
				if gene[2] in all_ORFs.keys():
					longest_ORF = ['', 0]
					for ORF in all_ORFs[gene[2]].keys():
						if len(all_ORFs[gene[2]][ORF][1]) > longest_ORF[1]:
							longest_ORF = [ORF, len(all_ORFs[gene[2]][ORF][1])]
					accession_gene_ORF[analysis[0]][gene[2]] = all_ORFs[gene[2]][longest_ORF[0]]
				else:
					short_BUSCO.append(gene[0])
					print gene[0], gene[2]
				
	BUSCO_genes_for_analysis = list(sets.Set(BUSCO_genes_for_analysis) - sets.Set(short_BUSCO))

	# make output directory
	process_name = "mkdir %s" % options.output
	process = subprocess.Popen(process_name, shell = True)
	process.wait()

	# export FASTA files for alignment
	PRANK_script_file = open(options.output + '/' + 'prank.sh', 'w')
	PRANK_script_file.write('#! /bin/bash' + '\n')
	PRANK_script_file.write('\n')

	for BUSCO_gene in BUSCO_genes_for_analysis:
		PRANK_script_file.write('prank -d=' + BUSCO_gene + '.fa -o=' + BUSCO_gene + '.PRANK.phy -f=phylips -DNA -codon' + '\n')
		FASTA_file = open(options.output + '/' + BUSCO_gene + '.fa', 'w')
		for species in BUSCO_gene_species[BUSCO_gene]:
			if BUSCO_gene in accession_BUSCO_transcript[species].keys():
				if accession_BUSCO_transcript[species][BUSCO_gene] in accession_gene_ORF[species].keys():
					FASTA_file.write('>' + accession_BUSCO_transcript[species][BUSCO_gene] + '\n')
					FASTA_file.write(accession_gene_ORF[species][accession_BUSCO_transcript[species][BUSCO_gene]][1] + '\n')
		FASTA_file.close()

	PRANK_script_file.close()

	return
	
if __name__ == '__main__':
	main()

