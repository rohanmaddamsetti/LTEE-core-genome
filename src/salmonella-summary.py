#!/usr/bin/python

## salmonella-summary.py by Rohan Maddamsetti

## This script prints out a csv table linking protein IDs to salmonella
## locus IDs.

## Only genes with consistent OMA and Hatcher panortholog calls are used.

## Usage: python salmonella-summary.py > ../results/salmonella-orthology-summary.csv 

from Bio import SeqIO

def main():
    
    salmonella_gene_hash = {}

    for seq_record in SeqIO.parse("/Users/Rohandinho/Desktop/Projects/LTEE-core-genome/data/salmonella.gbk", "genbank"):
        for feat in seq_record.features:
            if feat.type == "CDS":
                cur_locus_tag = feat.qualifiers['locus_tag'][0]
                cur_protein_id = feat.qualifiers['protein_id'][0]
                salmonella_gene_hash[cur_locus_tag] = cur_protein_id

    ecoli_gene_hash = {}
    core_summary = open("/Users/Rohandinho/Desktop/Projects/LTEE-core-genome/results/core-summary.csv", "r")
    for line in core_summary:
        data = line.split(",")
        prot_id = data[0]
        locus_tag = data[1]
        ecoli_gene_hash[locus_tag] = prot_id

    oma_ortholog_hash = {}
    oma_orthology = open("/Users/Rohandinho/Desktop/Projects/LTEE-core-genome/data/OMA-REL606-LT2-orthologs.txt", "r")
    for line in oma_orthology:
        data = line.split("\t")
        if data[2] == "1:1":  # 1:1 orthology between Salmonella and REL606
            ecoli_tag = data[0]
            salmonella_tag = data[1]
            oma_ortholog_hash[ecoli_tag] = salmonella_tag

    hatcher_ortholog_hash = {}
    hatcher_orthology = open("/Users/Rohandinho/Desktop/Projects/LTEE-core-genome/data/NC_012967-Salmonella-RBB--evalue.rbb-panorthologs.txt", "r")
    for line in hatcher_orthology:
        data = line.split("\t")
        ecoli_prot_tag = data[0][22:] # lazy way to slice the strings.
        salmonella_tag = data[1][21:32]
        hatcher_ortholog_hash[ecoli_prot_tag] = salmonella_tag

    ## only use genes where OMA and Hatcher panortholog calls are identical.
    header = ["locus_tag", "Protein_id", "Salmonella_locus_tag", "Salmonella_Protein_id"]
    print ",".join(header)

    #count = 0
    for k in sorted(ecoli_gene_hash.keys()):
        if oma_ortholog_hash.has_key(k):
            v1 = oma_ortholog_hash[k]
        else:
            continue
        if hatcher_ortholog_hash.has_key(ecoli_gene_hash[k]):
            v2 = hatcher_ortholog_hash[ecoli_gene_hash[k]]
        else:
            continue
        if not salmonella_gene_hash.has_key(v1):
            continue
        if salmonella_gene_hash[v1] != v2:
            continue
        #count += 1
        row = [k, ecoli_gene_hash[k], v1, v2] 
        print ",".join(str(v) for v in row)
    #print count
if __name__ == "__main__":
    main()
