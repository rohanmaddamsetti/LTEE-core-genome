#!/usr/bin/python

## tabulateMutations.py by Rohan Maddamsetti.
## This script takes a set of annotated genomediff files,
## and it produces a csv file tabulating locus_tag, the gene, the genome source,
## if the genome is a mutator, and dN and dS for that locus.
## IMPORTANT: This script ignores the araA marker mutation so that it isn't spuriously analyzed as a
## nonsynonymous mutation that evolves in parallel.
## The araA marker is at position 70867, and is a T->C D92G mutation.

## Usage: python tabulateMutations.py --non-mutator ../data/annotated_non-mutator_40K_diffs/*.gd --mutator ../data/annotated_mutator_40K_diffs/*.gd

import argparse
from os.path import basename
from csv import writer
import pprint
import sys

def parse_annotated_gd(gdfile, genomedict, ismutator=True):
	'''This function updates genomedict with the data contained in gdfile.'''		
	if ismutator:
		mutator = "TRUE"
	else:
		mutator = "FALSE"
	##get the name of the genome.	
	genome_name = basename(gdfile).split('.')[0]
	gd_handle = open(gdfile)
	genomedict[genome_name] = {}
	for line in gd_handle:
                ## check if mutation is araA marker--if so, skip.
                if "gene_name=araA" in line and "aa_position=92" in line and "snp_type=nonsynonymous" in line:
                        continue
		if "snp_type=synonymous" in line or "snp_type=nonsynonymous" in line:
			data = line.split("\t")
                        gene_name = None
                        locus_tag = None
			for elt in data:
                                if elt.startswith("gene_name="):
                                        gene_name = elt.split('=')[1]
                                elif elt.startswith("locus_tag="):
                                        locus_tag = elt.split('=')[1]
			##save these data in some data structure.
			gene_dict = genomedict[genome_name]
			if locus_tag not in gene_dict:
				##initialize the entry (without the gene name entry yet)
				gene_dict[locus_tag] = dict.fromkeys(["dN", "dS"], 0)
			## Now update the entry.
	 		genedata = gene_dict[locus_tag]
	 		genedata["gene"] = gene_name
			genedata["mutator"] = mutator
			if "snp_type=synonymous" in line:
				genedata["dS"] = genedata["dS"] + 1
			elif "snp_type=nonsynonymous" in line:
				genedata["dN"] = genedata["dN"] + 1

def main(args):
	mutator_gds = args.mutator
	non_mutator_gds = args.non_mutator
	genome_dict = {}
        if mutator_gds:
                for gd_file in mutator_gds:
                        parse_annotated_gd(gd_file, genome_dict, ismutator=True)
	if non_mutator_gds:
                for gd_file in non_mutator_gds:
                        parse_annotated_gd(gd_file, genome_dict, ismutator=False)
	## Dump the data into a csv file.
	output = writer(open("../results/ltee_mutations.csv","wb"))
	header = ["locus_tag", "gene", "genome", "mutator", "dN", "dS"]
	output.writerow(header)
	for this_genome, this_gene_dict in genome_dict.iteritems():
	 	for this_locus, fields in this_gene_dict.iteritems():
	 		rowdata = [this_locus, fields["gene"], this_genome, fields["mutator"], 
					   fields["dN"], fields["dS"]]
			output.writerow(rowdata)

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("--non-mutator", help="list of non-mutator *.gd files", nargs='*', dest="non_mutator")
	parser.add_argument("--mutator", help="list of mutator *.gd files", nargs='*',dest="mutator")
	args = parser.parse_args()
	main(args)
