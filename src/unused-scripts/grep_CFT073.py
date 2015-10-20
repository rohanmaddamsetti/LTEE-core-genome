#!/usr/bin/python

# grep_CFT073.py by Rohan Maddamsetti

# This script pulls out the relevant E. coli CFT073 proteins under positive selection in the 
# Chattopadhyay 2012 dataset.

## NOTE: This uses an old (2002) Genbank file--parsing is for this file only.

# Usage: python grep_CFT073.py > ../results/relevant-CFT073-proteins.fasta

from Bio import SeqIO

def main():
    relevant_genes = []
    gene_file = open("../data/Chattopadhyay-core-gene-ids.txt")
    for i,line in enumerate(gene_file):
        line = line.strip()
        if i == 0:
            continue
        relevant_genes.append(line)
    #print relevant_genes

    CFT073_gb_file = open("../data/AE014075.gbk", "r")
    record =  SeqIO.parse(CFT073_gb_file,"genbank").next()
    for feat in record.features:
        if feat.type == "CDS":
            cur_gene = feat.qualifiers['gene'][0]
            prot_sequence = feat.qualifiers['translation'][0]
            if cur_gene in relevant_genes:
                print '>'+cur_gene+'=CFT073-gene |Escherichia coli CFT073|AE014075.gbk'
                print prot_sequence
                print
                
main()
