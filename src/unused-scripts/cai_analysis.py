#!/usr/bin/python

## cai_analysis.py by Rohan Maddamsetti.

## This script requires the EMBOSS programs cai and cusp to be in the system path.

## 1) Get all genes annotated as either "30S ribosomal protein" or "50S ribosomal protein"
##    in the REL606 genome.
## 2) Use these sequences as an input to cusp to make a reference codon usage table.
## 3) For all CDS in REL606, calculate CAI with the emboss program cai.
## 4) print out CSV formatted data (locus_tag, CAI) to stdout.

## Usage: python cai_analysis.py > ../data/cai.csv

from Bio import SeqIO
from Bio.Seq import Seq
from sys import exit
import subprocess
from os import mkdir, listdir
from os.path import isdir, basename
from string import split

def getRibosomalReferences(rel606_genbank_path, out_path):
	outfile = open(out_path, "w")
	for genome in SeqIO.parse(rel606_genbank_path, "genbank"):
		for feature in genome.features:
			if feature.type != "CDS":
				continue
			else:
				try:
					product = feature.qualifiers["product"][0]
					if product.startswith("30S ribosomal protein") or \
						product.startswith("50S ribosomal protein"):
						locus_tag = feature.qualifiers["locus_tag"][0]
						#print locus_tag
						start = feature.location.start.position
						end = feature.location.end.position
						annotated_translation = feature.qualifiers["translation"][0]
						sequence = genome.seq[start:end]
						try:
							translation = sequence.translate(table="Bacterial",
															 to_stop=True,cds=True)
						except:
							sequence = sequence.reverse_complement()
							translation = sequence.translate(table="Bacterial",
															 to_stop=True, cds=True)

						outfile.write(">" + locus_tag + "\n" + str(sequence) + "\n\n")
						#print sequence, "\n\n"
						#print translation, "\n", annotated_translation, "\n"
				except KeyError:
					continue
	outfile.close()

def run_cusp(reference_file_path, outfile_path):
	## the command is this: cusp -sequence ribosomal_references.fasta -outfile ribosomal_references.cusp
	args = ["cusp", "-sequence", reference_file_path, "-outfile", outfile_path]
	p = subprocess.Popen(args)

def print_cds_sequences(relpath, output_dir):
	for genome in SeqIO.parse(relpath, "genbank"):
		for feature in genome.features:
			if feature.type != "CDS":
				continue
			else:
				try:
					locus_tag = feature.qualifiers["locus_tag"][0]
					#print locus_tag
					annotated_translation = feature.qualifiers["translation"][0]
					start = feature.location.start.position
					end = feature.location.end.position
					sequence = genome.seq[start:end]
					try:
						translation = sequence.translate(table="Bacterial",
														 to_stop=True,cds=True)
					except:
						sequence = sequence.reverse_complement()
						try:
							translation = sequence.translate(table="Bacterial",
														 to_stop=True, cds=True)
						except: #This checking fails on CDS containing selenocysteine.
							continue
					outfile = open(output_dir+locus_tag+".fasta", "w")
					outfile.write(">" + locus_tag + "\n" + str(sequence) + "\n\n")
					outfile.close()
				except KeyError:
					continue
    
def run_cai(input_sequence_path, table_file_path, output_cai_path):
	##print input_sequence_path
	outfile_path = split(basename(input_sequence_path), ".")[0] + ".cai"
	##print outfile_path
	args = ["cai", "-seqall", input_sequence_path, "-cfile", table_file_path, "-outfile", output_cai_path+outfile_path]
	##print args
	p = subprocess.Popen(args)
    
def main():
	rel606path = "../data/REL606.gbk"
	reference_path = "../data/ribosomal_references.fasta"
	getRibosomalReferences(rel606path, reference_path)
	## Now, run cusp on the ribosomal references.
	cusp_table_path = "../data/ribosomal_references.cusp"
	run_cusp(reference_path, cusp_table_path)
	## make a directory containing every CDS sequence in REL606.gbk.
	cds_path = "../data/cds_sequences/"
	if not isdir(cds_path):
		mkdir(cds_path)
	print_cds_sequences(rel606path, cds_path)
    ## Now, run cai on every CDS sequence in REL606.
	cai_path = "../data/cai_values/"
	if not isdir(cai_path):
		mkdir(cai_path)
	for fasta in listdir(cds_path):
		run_cai(cds_path+fasta, cusp_table_path, cai_path)
	## Now, make a csv file with locus_tag, CAI values.
	print "locus_tag, CAI"
	for cai_file in listdir(cai_path):
		cai_f = open(cai_path+cai_file, "r")
		for line in cai_f:
			data = line.split()
			locus_tag = data[1]
			cai = data[3]
			print locus_tag + ",", cai
		cai_f.close()

if __name__ == "__main__":			
	main()
