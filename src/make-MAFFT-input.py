#!/usr/bin/python

# make-MAFFT-input.py by Rohan Maddamsetti.

# This script generates FASTA files that MAFFT-linsi will use
# to produce protein alignments.

import re
from os import listdir

def make_FASTAdict(path):
    sequence_file = open(path, "r")
    header = ""
    fasta_seq = ""
    protein_id = ""
    FASTAdict = {}
    for line in sequence_file:
        if line.startswith(">"):
            if len(fasta_seq) > 0: # dump old entry into the dict.
                FASTAdict[protein_id] = (header, fasta_seq)
            header = line
            fasta_seq = ""
            protein_id = ""
            protein_id = re.search("\[protein_id=(.+?)\]", header).group(1)
        else:
            fasta_seq = fasta_seq + line
    ## at the end of the loop, dump last entry into the dict.
    if len(fasta_seq) > 0: # dump old entry into the dict.
        FASTAdict[protein_id] = (header, fasta_seq)
    return(FASTAdict)

def make_salmonella_ecoli_MAFFT_inputs():
    localdir = "/Users/Rohandinho/Desktop/Projects/LTEE-core-genome/"
    outputdir = "data/MAFFT-input/"
    orthology_summary = open(localdir+"results/salmonella-orthology-summary.csv","r")

    REL606_prot_path = localdir + "data/CDS-sequences/NC_012967.txt"
    salmonella_prot_path = localdir + "data/CDS-sequences/salmonella-CDS-sequences.txt"

    REL606_prot_dict = make_FASTAdict(REL606_prot_path)
    salmonella_prot_dict = make_FASTAdict(salmonella_prot_path)

    for i, line in enumerate(orthology_summary):
        if i == 0:
            continue
        data = [x.strip() for x in line.split(',')]
        locus_tag, protein_id, salmonella_locus_tag, salmonella_protein_id = data

        ## make protein input FASTA file.
        prot_output_dir = localdir + outputdir + "protein/divergence/"
        prot_output = open(prot_output_dir + locus_tag + "_ecoli-salmonella_prot.txt", "w")

        prot_output.write(REL606_prot_dict[protein_id][0])
        prot_output.write(REL606_prot_dict[protein_id][1])
        prot_output.write(salmonella_prot_dict[salmonella_protein_id][0])
        prot_output.write(salmonella_prot_dict[salmonella_protein_id][1])
        prot_output.close()

def make_panortholog_MAFFT_inputs():
    localdir = "/Users/Rohandinho/Desktop/Projects/LTEE-core-genome/"
    outputdir = "data/MAFFT-input/"
    prot_path = "data/CDS-sequences/"

    ## make a dictionary of REL606 protein_ids to locus_tags.
    REL606dict = {}
    orthology_summary = open(localdir+"results/core-summary.csv","r")
    #orthology_summary = open(localdir+"results/15-genomes-results/core-summary.csv","r")
    for i, line in enumerate(orthology_summary):
        data = [x.strip() for x in line.split(',')]
        protein_id = data[0]
        locus_tag = data[1]
        REL606dict[protein_id] = locus_tag

    ## make dictionaries of all FASTAdicts.
    all_prot_FASTAdicts = {}
    full_prot_path = localdir + prot_path
    for seqfile in listdir(full_prot_path):
        if seqfile.startswith("NC"): #skip salmonella sequences.
            gbk_id = seqfile.split(".")[0]
            all_prot_FASTAdicts[gbk_id] = make_FASTAdict('/'.join([full_prot_path,seqfile]))

    orthology_file = open(localdir+"data/60-genomes-hatcher-results/pan-9-59.tmp.txt")
    #orthology_file = open(localdir+"data/15-genomes-hatcher-results/pan-1-14.tmp.txt")
    for i, line in enumerate(orthology_file):
        data = [x.strip() for x in line.split()]
        REL606protein_id =  re.search("NC_012967$lcl|NC_012967\.1_cdsid_(.+?)\t", line).group(1)
        prot_output_dir = localdir + outputdir + "protein/polymorphism/"
        #prot_output_dir = localdir + outputdir + "protein/15-genomes-polymorphism/"
        prot_output = open(prot_output_dir + REL606dict[REL606protein_id] + "_ecoli-orthologs.txt", "w")
        for elt in data:
            m2 = re.search("(.+?)\$lcl\|(.+?)_cdsid_(.+)", elt)
            gbk = m2.group(1)
            protein_id = m2.group(3)
            prot_output.write(all_prot_FASTAdicts[gbk][protein_id][0])
            prot_output.write(all_prot_FASTAdicts[gbk][protein_id][1])
        prot_output.close()

def make_tree_MAFFT_inputs():
    localdir = "/Users/Rohandinho/Desktop/Projects/LTEE-core-genome/"
    outputdir = "data/tree-MAFFT-input/"
    prot_path = "data/CDS-sequences/"
    tree_orthology_summary = open(localdir+"results/tree/tree-orthology-summary.csv","r")

    salmonella_prot_path = localdir + "data/CDS-sequences/salmonella-CDS-sequences.txt"
    salmonella_prot_dict = make_FASTAdict(salmonella_prot_path)
    ## make dictionaries of all FASTAdicts for E. coli strains.
    strain_prot_FASTAdicts = {}
    full_prot_path = localdir + prot_path
    for seqfile in listdir(full_prot_path):
        if seqfile.startswith("NC"): #skip salmonella sequences, handle separately.
            gbk_id = seqfile.split(".")[0]
            strain_prot_FASTAdicts[gbk_id] = make_FASTAdict('/'.join([full_prot_path,seqfile]))

    ## dict of REL606 protein_ids to (locus_tag, salmonella_locus_tag,salmonella_protein_id)
    ## for tree genes of interest.
    REL606_treedict = {}
    for i, line in enumerate(tree_orthology_summary):
        if i == 0:
            continue
        data = [x.strip() for x in line.split(',')]
        locus_tag, protein_id, salmonella_locus_tag, salmonella_protein_id, note, start, length, panortholog = data
        REL606_treedict[protein_id] = (locus_tag, salmonella_locus_tag, salmonella_protein_id)

    panorthology_file = open(localdir+"data/60-genomes-hatcher-results/pan-9-59.tmp.txt")
    for i,line in enumerate(panorthology_file):
        data = [x.strip() for x in line.split()]
        REL606protein_id =  re.search("NC_012967$lcl|NC_012967\.1_cdsid_(.+?)\t", line).group(1)
        ## skip genes that aren't of interest for tree building.
        if REL606protein_id not in REL606_treedict:
            continue
        locus_tag, salmonella_locus_tag, salmonella_protein_id = REL606_treedict[REL606protein_id]
        prot_output_dir = localdir + outputdir
        prot_output = open(prot_output_dir + locus_tag + ".txt", "w")
        ## write Salmonella ortholog into the MAFFT input file.
        salmonella_header,salmonella_seq = salmonella_prot_dict[salmonella_protein_id]
        prot_output.write(salmonella_header)
        prot_output.write(salmonella_seq)
        ## write E. coli strain panorthologs into the MAFFT input file.
        for elt in data:
            m2 = re.search("(.+?)\$lcl\|(.+?)_cdsid_(.+)", elt)
            gbk = m2.group(1)
            protein_id = m2.group(3)
            strain_header, strain_seq = strain_prot_FASTAdicts[gbk][protein_id]
            prot_output.write(strain_header)
            prot_output.write(strain_seq)
        prot_output.close()



def main():
    #make_salmonella_ecoli_MAFFT_inputs()
    #make_panortholog_MAFFT_inputs()
    make_tree_MAFFT_inputs()


main()
