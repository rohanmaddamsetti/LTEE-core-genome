#!/usr/bin/python

## protein-site-diversity.py by Rohan Maddamsetti.

## Usage: python protein-site-diversity.py > ../results/protein-site-diversity.csv

## OR python protein-site-diversity.py > ../results/protein-site-divergence.csv
## when calculating from Salmonella data. Be sure to check main before running to
## comment out/run the proper function calls.

# Make a list of dN in non-mutators:
# [gene, codon, genome]

## IMPORTANT: This script ignores the araA marker mutation so that it isn't spuriously analyzed## as a  nonsynonymous mutation that evolves in parallel.
## The araA marker is at position 70867, and is a T->C D92G mutation.

## IMPORTANT TOO: This script ignores a recD secondary mutation found in Ara+ lineages
## (this mutation occurs in REL607 and its descendants).


from __future__ import division
from os import listdir, rename, getcwd, remove
import os.path
import re
import subprocess
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from pprint import pprint

def make_dN_dict(diff_directory):
    dN_dict = {}
    for diff in [x for x in listdir(diff_directory) if x.endswith('.gd')]:
        cur_path = os.path.join(diff_directory, diff)
        cur_handle = open(cur_path)
        genome = diff.partition(".")[0]
        nonsynonymous_count = 0
        for line in cur_handle:
            if line.startswith("SNP"):
                ## check if mutation is araA marker--if so, skip.
                if "gene_name=araA" in line and "aa_position=92" in line and "snp_type=nonsynonymous" in line:
                    continue
                ## check if mutation is recD secondary mutation in 607--if so, skip.
                if "gene_name=recD" in line and "aa_position=10" in line and "snp_type=nonsynonymous" in line:
                    continue
                if "snp_type=nonsynonymous" in line:
                    nonsynonymous_count = nonsynonymous_count + 1
                    fields = line.split("\t")
                    #print fields
                    for i in fields:
                        if i.startswith("gene_name"):
                            gene_name = i.split("=")[1]
                        if i.startswith("locus_tag"):
                            locus_tag = i.split("=")[1]
                        if i.startswith("aa_position"):
                            aa_position = int(i.split("=")[1])
                    if locus_tag in dN_dict:
                        dN_dict[locus_tag]["LTEE_dN"] += 1
                        dN_dict[locus_tag]["aa_position"] |= {aa_position}
                    else:
                        dN_dict[locus_tag] = {
                            "locus_tag":locus_tag,
                            "gene_name":gene_name,
                            "aa_position":{aa_position},
                            "LTEE_dN":1
                        }
    return dN_dict

def calculate_protein_site_diversity(dN_dict,diverge=False):
    '''Calculate Nei's nucleotide diversity, modified for peptide sequences. This function
    runs through proteins that evolved in the LTEE, and calculates diversity at those
    amino acid positions that mutated and diversity at the rest of the protein.
    Are the sites that mutated in the LTEE more or less variable than the rest of the protein
    in nature? '''
    folder1="protein"
    if diverge == False:
        pattern = '^(.+)_ecoli-orthologs'
        folder2 = "15-genomes-polymorphism"
    elif diverge == True:
        pattern = '^(.+)_ecoli-salmonella_prot'
        folder2 = "divergence"

    localdir = "/Users/Rohandinho/Desktop/Projects/LTEE-core-genome/"
    alignments_dir = os.path.join(localdir, "results/MAFFT-output/", folder1, folder2)
    print("locus_tag, mutated_site_diversity, rest_protein_diversity, gene_name, LTEE_dN")
    for in_aln in listdir(alignments_dir):
        if in_aln.startswith('ECB'):
            filename = in_aln.partition('.')[0]
            locus_tag = re.search(pattern, filename).group(1)
            if locus_tag not in dN_dict:
                continue
            else:
                alignment = AlignIO.read(os.path.join(alignments_dir,in_aln), "fasta")
                n = len(alignment)
                mutated_site_diversity = float(0)
                rest_protein_diversity = float(0)
                d_denominator = n*(n-1)/2 # number of comparisons is n choose 2.
                mutated_sites = len(dN_dict[locus_tag]["aa_position"])
                rest_sites  = alignment.get_alignment_length() - mutated_sites
                for j in range(n):
                    for i in range(j):
                        sub_aln = alignment[i:j+1:j-i]
                        mutated_site_differences = 0
                        rest_site_differences = 0
                        for x in range(alignment.get_alignment_length()):
                            cur_site = sub_aln[:,x]
                            if len(set(cur_site)) > 1:
                                if (x+1) in dN_dict[locus_tag]["aa_position"]:
                                    mutated_site_differences += 1
                                else:
                                    rest_site_differences += 1
                        pair_mutated_site_diversity = mutated_site_differences/mutated_sites
                        pair_rest_site_diversity = rest_site_differences/rest_sites
                        mutated_site_diversity += pair_mutated_site_diversity
                        rest_protein_diversity += pair_rest_site_diversity
                d_denominator = n*(n-1)/2 # number of comparisons is n choose 2.
                mutated_site_diversity = mutated_site_diversity/d_denominator
                rest_protein_diversity = rest_protein_diversity/d_denominator
                print(','.join([locus_tag, str(mutated_site_diversity), str(rest_protein_diversity), dN_dict[locus_tag]["gene_name"], str(dN_dict[locus_tag]["LTEE_dN"])]))


def main():
    proj_dir = "/Users/Rohandinho/Desktop/Projects/LTEE-core-genome/"
    diff_dir = proj_dir + "data/annotated_non-mutator_50K_diffs"
    dN_dict = make_dN_dict(diff_dir)
    #calculate_protein_site_diversity(dN_dict,diverge=False)
    calculate_protein_site_diversity(dN_dict,diverge=True)


main()
