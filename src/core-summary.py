#!/usr/bin/python

## core-summary.py by Rohan Maddamsetti

## This script prints out a csv table that classifies genes in REL606
## as being either panorthologs or not, based on the analysis that
## Phil Hatcher ran for me.

## Unfortunately, this script uses two different reference genbank files:
## the REL606 found in genbank (NCBI-REL606.gbk), and the reference
## used to generate the newest set of genomediffs from the LTEE
## (REL606.6.gbk).

## Usage: python core-summary.py > ../results/core-summary.csv 

from Bio import SeqIO

def main():
    
    core_genes = []
    f1 = open("../data/15-genomes-hatcher-results/remainingGenes.15.txt")
    #f1 = open("../data/60-genomes-hatcher-results/remainingGenes.txt")
    for line in f1:
        core_genes.append(line.strip())
    f1.close()

    fpath = "/Users/Rohandinho/Desktop/Projects/LTEE-core-genome/data/"
    locus_hash = {}
    ## Get the 'note' data from REL606.6.gbk.
    for seq_record in SeqIO.parse(fpath+"REL606.6.gbk", "genbank"):
        for feat in seq_record.features:
            if feat.type == "CDS":
                cur_locus_tag = feat.qualifiers['locus_tag'][0]
                try:
                    notes = feat.qualifiers['note'][0].split(";")
                    for i in notes:
                        if i.startswith('b'):
                            cur_note = i
                except KeyError:
                    cur_note = "NA"
                locus_hash[cur_locus_tag] = cur_note

    gene_hash = {}
    gene_total = 0 # sanity check


    ## Tabulate length, protein_id, and locus_tag for all genes in REL606
    ## reference genbank file.
    for seq_record in SeqIO.parse(fpath+"NCBI-REL606.gbk",
                                  "genbank"):
        for feat in seq_record.features:
            if feat.type == "CDS":
                #print feat.qualifiers
                gene_total = gene_total + 1
                cur_protein_id = feat.qualifiers['protein_id'][0]
                cur_locus_tag = feat.qualifiers['locus_tag'][0]
                cur_start_pos = str(feat.location.start+1)
                cur_length = abs(int(feat.location.start) - int(feat.location.end))
                try:
                    cur_note = locus_hash[cur_locus_tag]
                except KeyError:
                    continue
                if cur_protein_id in core_genes:
                    gene_hash[cur_protein_id] = (cur_locus_tag, cur_note, cur_start_pos, cur_length,"TRUE")
                else:
                    gene_hash[cur_protein_id] = (cur_locus_tag, cur_note, cur_start_pos, cur_length,"FALSE")

    header = ["Protein_id", "locus_tag", "note", "Start position", "length", "panortholog"]
    print(",".join(header))
    for k in sorted(gene_hash.keys()):
        row = [k] + list(gene_hash[k])
        print(",".join(str(v) for v in row))

if __name__ == "__main__":
    main()
