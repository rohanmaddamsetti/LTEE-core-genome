#!/usr/bin/python

'''calculate-diversity.py by Rohan Maddamsetti

USAGE:
See main function for usage. Pretty ugly right now involving commenting out lines.

Align orthologous genes between REL606 and S. typhimurium LT2.
Count the number of sites that have diverged between REL606 and Salmonella
for each ortholog.

Align orthologous genes among the 60 sampled E. coli.
Count the number of polymorphic sites in each panortholog.

In both cases, normalize by gene length and the number of pairwise comparisons made to
calculate sequence diversity.

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
import operator
import itertools

# compute hamming distance of two sequences
def hamming(s1,s2):
    return sum(map(operator.ne,s1,s2))

def calculate_diversity(diverge=True):
    '''Calculate Nei's nucleotide diversity, modified for peptide sequences. '''

    folder1 = "protein"
    if diverge == True:
        pattern = '^(.+)_ecoli-salmonella'
        folder2 = "divergence"
    else:
        pattern = '^(.+)_ecoli-orthologs'
        folder2 = "15-genomes-polymorphism"

    localdir = "/Users/Rohandinho/Desktop/Projects/LTEE-core-genome/"
    alignments_dir = os.path.join(localdir, "results/MAFFT-output/", folder1, folder2)
    print("locus_tag, sites, diversity")
    for in_aln in listdir(alignments_dir):
        if in_aln.startswith('ECB'):
            filename = in_aln.partition('.')[0]
            locus_tag = re.search(pattern, filename).group(1)
            try:
                alignment = AlignIO.read(os.path.join(alignments_dir,in_aln), "fasta")
                n = len(alignment)
                diversity = float(0)
                d_denominator = n*(n-1)/2 # number of comparisons is n choose 2.
                sites = alignment.get_alignment_length()
                for i,j in itertools.combinations(range(n),2):
                    sub_aln = alignment[i:j+1:j-i]
                    top = sub_aln[0]
                    bottom = sub_aln[-1]
                    print(top.seq)
                    print(bottom.seq)
                    differences = hamming(top.seq,bottom.seq)
                    pair_diversity = differences/sites
                    diversity += pair_diversity
                diversity = diversity/d_denominator
                print(','.join([locus_tag, str(sites), str(diversity)]))
            except ValueError: ## if MAFFT failed
                continue

def main():

# DO THIS:    python calculate-diversity.py > ../results/ecoli-protein-diversity.csv
    calculate_diversity(diverge=False)

# DO THIS:    python calculate-diversity.py > ../results/ecoli-salmonella-protein-diversity.csv
    #calculate_diversity(diverge=True)

if __name__ == '__main__':
    main()
