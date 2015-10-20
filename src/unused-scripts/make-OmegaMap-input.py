#!/usr/bin/python

## This script copies all the alignment files in 
## results/MAFFT-output/nucleotide/polymorphism
## for genes that mutated in the LTEE at 40K
## to data/omegaMap-input/nucleotide/polymorphism.
## after stripping stop codons for omegaMap to work,
## and making sure that the sequence is mod 3 (a coding sequence).

from os import listdir
import os.path
from Bio import AlignIO
import sys
import re

## helper function to remove duplicates from a list.
def f7(seq):
    seen = set()
    seen_add = seen.add
    return [x for x in seq if x not in seen and not seen_add(x)]

def copy_alns(divergence=True):
    local_dir = "/Users/Rohandinho/Desktop/Projects/LTEE-core-genome/"
    
    folder = None
    if divergence:
        folder = "divergence"
    else:
        folder = "polymorphism"

    data_path= os.path.join(local_dir,"results/MAFFT-PRANK-output/nucleotide/", folder)
## the following is where the files will be put.
    omega_input_dir = os.path.join(local_dir, "data/omegaMap-input/", folder)

    # get the locus_tag of genes with mutations in 40K nonmutators.
    relevant_gene_file = open(os.path.join(local_dir,"results/ltee_mutations_40K.csv") , "r")
    relevant_gene_list = []
    for num, line in enumerate(relevant_gene_file):
        if num == 0:
            continue
        x =  line.split(",")
        if x[3] == "TRUE": # is a mutator
            continue
        elif x[3] == "FALSE": # is a non-mutator
            locus_tag = x[0]
            if not (locus_tag in relevant_gene_list):
                relevant_gene_list.append(locus_tag)
    relevant_gene_list = f7(relevant_gene_list)
    relevant_gene_list.sort()
    #print relevant_gene_list
    #print len(relevant_gene_list)

    for f in listdir(data_path):
        if f.startswith('ECB') and '.aln' in f:
            ## see if the gene mutated in the LTEE for polymorphism analysis.
            tag = re.search("(^ECB_[0-9]{5})_",f).group(1)
            if (folder == "divergence" ) or (tag in relevant_gene_list):
                filename = os.path.join(data_path,f)
                alignment = AlignIO.read(filename, "fasta")
                alignment_length = alignment.get_alignment_length()
                if alignment_length % 3 != 0:
                    print "Warning! alignment", f, "is not in frame!"
                    continue
                else:
                    # strip the terminal stop codon (assume that it's there).
                    alignment2 = alignment[:, :-3]
                    # write to file.
                    outfname = os.path.join(omega_input_dir, f) 
                    #print alignment2
                    #print outfname
                    outfile = open(outfname , "w")
                    AlignIO.write(alignment2, outfile, "fasta")
                    outfile.close()
            else:
                continue
        else:
            continue


def main():
    copy_alns(divergence=True)
    copy_alns(divergence=False)
main()
