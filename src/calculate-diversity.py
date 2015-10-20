#!/usr/bin/python

'''calculate-diversity.py by Rohan Maddamsetti

Usage: 
     python calculate-diversity.py > ../results/ecoli-salmonella-protein-diversity.csv
     python calculate-diversity.py > ../results/ecoli-protein-diversity.csv

Align orthologous genes between REL606 and S. typhimurium LT2.
Count the number of sites that have diverged between REL606 and Salmonella
for each ortholog.

Align orthologous genes among the 60 sampled E. coli.
Count the number of polymorphic sites in each panortholog.

'''
from __future__ import division
from math import log
from os import listdir, rename, getcwd, remove
import os.path
import re
import subprocess
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from pprint import pprint

def calculate_diversity(nucleotide=False, diverge=True):
    '''Calculate Nei's nucleotide diversity, modified for peptide sequences. '''
    if nucleotide == True:
        folder1="nucleotide"
    else:
        folder1="protein"
    if diverge == True:
        pattern = '^(.+)_ecoli-salmonella'
        folder2 = "divergence"
    else:
        pattern = '^(.+)_ecoli-orthologs'
        folder2 = "polymorphism"

    localdir = "/Users/Rohandinho/Desktop/Projects/LTEE-core-genome/"    
    alignments_dir = os.path.join(localdir, "results/MAFFT-output/", folder1, folder2)
    print "locus_tag, sites, diversity"
    for in_aln in listdir(alignments_dir):
        if in_aln.startswith('ECB'):
            filename = in_aln.partition('.')[0]
            locus_tag = re.search(pattern, filename).group(1)
            try:
                alignment = AlignIO.read(os.path.join(alignments_dir,in_aln), "fasta")
                n = len(alignment)
                #print alignment
                diversity = float(0)
                d_denominator = n*(n-1)/2 # number of comparisons is n choose 2.
                sites = alignment.get_alignment_length()
                for j in range(n):
                    for i in range(j):
                        sub_aln = alignment[i:j+1:j-i]
                        differences = 0
                        for x in range(sites):
                            cur_site = sub_aln[:,x]
                            if len(set(cur_site)) > 1:
                                differences = differences + 1
                        pair_diversity = differences/sites
                        diversity += pair_diversity
                d_denominator = n*(n-1)/2 # number of comparisons is n choose 2.
                #print n
                #print d_denominator
                #print "comparisons: ", d_denominator
                diversity = diversity/d_denominator
                #print "diversity: ", diversity
                print ','.join([locus_tag, str(sites), str(diversity)])
            except ValueError: ## if MAFFT failed
                continue

def calculate_differences(nucleotide=False, diverge=True):
    if nucleotide == True:
        folder1="nucleotide"
    else:
        folder1="protein"
    if diverge == True:
        pattern = '^(.+)_ecoli-salmonella'
        folder2 = "divergence"
    else:
        pattern = '^(.+)_ecoli-orthologs'
        folder2 = "polymorphism"
    localdir = "/Users/Rohandinho/Desktop/Projects/LTEE-core-genome/"    
    alignments_dir = os.path.join(localdir, "results/MAFFT-output/", folder1, folder2)
    print "locus_tag, sites, differences"
    for in_aln in listdir(alignments_dir):
        if in_aln.startswith('ECB'):
            filename = in_aln.partition('.')[0]
            locus_tag = re.search(pattern, filename).group(1)
            try:
                alignment = AlignIO.read(os.path.join(alignments_dir,in_aln), "fasta")
                sites = alignment.get_alignment_length()
                differences = 0
                for i in range(sites):
                    cur_site = alignment[:,i]
                    #print cur_site
                    if len(set(cur_site)) > 1:
                        differences = differences + 1
                print ','.join([locus_tag, str(sites), str(differences)])
            except ValueError: ## if MAFFT failed
                continue
def main():
    # print out csv file with locus_tag, sites, diversity.
    calculate_diversity(nucleotide=False, diverge=False)
    #calculate_diversity(nucleotide=False, diverge=True)    

if __name__ == '__main__':
    main()
