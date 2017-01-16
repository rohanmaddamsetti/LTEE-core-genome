#!/usr/env python

## filter-membrane-proteins.py by Rohan Maddamsetti.

## this script is for preparing data for the solvent accessibility
## analysis.

## First, get all REL606 genes with observed dN
## (in non-mutators or hypermutators).
## Loop through the UniProt annotation,
## and mark these genes with no membrane annotation (0),
## genes with membrane annotation (1)
## and genes not found in the UniProt annotation (2).
## End when all genes have been annotated or when
## the whole UniProt file has been read.

## NOTE: ONLY proteins with LTEE hits in Tenaillon et al. 2016 are annotated!

## Usage: python filter-membrane-proteins.py > ../results/membrane-annotation.csv

import pandas as pd

def main():

    tdata2 = pd.read_csv("../data/Tenaillon-data/nature18959-s2.csv")
    tdata3 = pd.read_csv("../data/Tenaillon-data/nature18959-s3.csv")

    genes2 = tdata2.loc[tdata2['Observed nonsynonymous mutation'] != 0]['Gene name']
    genes3 = tdata3.loc[tdata3['Observed nonsynonymous mutation'] != 0]['Gene name']
    relevant_genes = set(genes2) | set(genes3)

    gene_order = {x:y for x,y in zip(tdata2['Gene order'],tdata2['Gene name'])}

    gene_dict = {x:2 for x in relevant_genes}
    membrane_annotation = open("../data/Ecoli-uniprot.txt")
    for line in membrane_annotation:
        if len(relevant_genes) == 0:
            break
        line = line.strip()
        dat = line.split()
        for x in relevant_genes:
            if x in line:
                if 'membrane' in line and 'Cytoplasm' not in line:
                    gene_dict[x] = 1
                elif 'transport' in line or 'Transport' in line:
                    gene_dict[x] = 1
                else:
                    gene_dict[x] = 0
                relevant_genes.remove(x)
                break

    relevant_genes = gene_dict.keys()
    print('Gene order ,Gene name, membrane annotation')
    for y in sorted(gene_order.keys()):
        z = gene_order[y]
        if z in relevant_genes:
            print(str(y)+','+str(z)+','+str(gene_dict[z]))

main()
