#!/usr/bin/python

# make-MAFFT-input.py by Rohan Maddamsetti.

# This script generates FASTA files that MAFFT-linsi will use
# to produce protein alignments.
# lines relevant to making nucleotide alignments have been commented out.


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

    #REL606_nuc_path = localdir + "data/nucleotide-sequences/NC_012967.txt"
    REL606_prot_path = localdir + "data/CDS-sequences/NC_012967.txt"
    #salmonella_nuc_path = localdir + "data/nucleotide-sequences/salmonella-nucleotide.txt"
    salmonella_prot_path = localdir + "data/CDS-sequences/salmonella-CDS-sequences.txt"

    #REL606_nuc_dict = make_FASTAdict(REL606_nuc_path)
#    for k in sorted(REL606_nuc_dict.keys()):
#        print k, REL606_nuc_dict[k]
    REL606_prot_dict = make_FASTAdict(REL606_prot_path)
    #salmonella_nuc_dict = make_FASTAdict(salmonella_nuc_path)
    salmonella_prot_dict = make_FASTAdict(salmonella_prot_path)

    for i, line in enumerate(orthology_summary):
        if i == 0:
            continue
        data = [x.strip() for x in line.split(',')]
        locus_tag, protein_id, salmonella_locus_tag, salmonella_protein_id = data
        
        ## make nucleotide input FASTA file.
        #nuc_output_dir = localdir + outputdir + "nucleotide/divergence/"
        #nuc_output = open(nuc_output_dir + locus_tag + "_ecoli-salmonella_nuc.txt", "w")
        #print REL606_nuc_dict[protein_id][0]
        #print REL606_nuc_dict[protein_id][1]
        #print salmonella_nuc_dict[salmonella_protein_id][0]
        #print salmonella_nuc_dict[salmonella_protein_id][1]
   
        #nuc_output.write(REL606_nuc_dict[protein_id][0])
        #nuc_output.write(REL606_nuc_dict[protein_id][1])
        #nuc_output.write(salmonella_nuc_dict[salmonella_protein_id][0])
        #nuc_output.write(salmonella_nuc_dict[salmonella_protein_id][1])
        #nuc_output.close()

        ## make protein input FASTA file.
        prot_output_dir = localdir + outputdir + "protein/divergence/"
        prot_output = open(prot_output_dir + locus_tag + "_ecoli-salmonella_prot.txt", "w")

        #print REL606_prot_dict[protein_id][0]
        #print REL606_prot_dict[protein_id][1]
        #print salmonella_prot_dict[salmonella_protein_id][0]
        #print salmonella_prot_dict[salmonella_protein_id][1]

        prot_output.write(REL606_prot_dict[protein_id][0])
        prot_output.write(REL606_prot_dict[protein_id][1])
        prot_output.write(salmonella_prot_dict[salmonella_protein_id][0])
        prot_output.write(salmonella_prot_dict[salmonella_protein_id][1])
        prot_output.close()
    
def make_panortholog_MAFFT_inputs():
   localdir = "/Users/Rohandinho/Desktop/Projects/LTEE-core-genome/"
   outputdir = "data/MAFFT-input/"
   #nuc_path = "data/nucleotide-sequences/"
   prot_path = "data/CDS-sequences/"

   # make a dictionary of REL606 protein_ids to locus_tags.
   REL606dict = {}
   orthology_summary = open(localdir+"results/core-summary.csv","r")
   for i, line in enumerate(orthology_summary):
       if i == 0:
           continue
       data = [x.strip() for x in line.split(',')]
       protein_id = data[0]
       locus_tag = data[1]
       REL606dict[protein_id] = locus_tag

   # make dictionaries of all FASTAdicts.
   #all_nuc_FASTAdicts = {}   
   #full_nuc_path = localdir + nuc_path
   #for seqfile in listdir(full_nuc_path):
   #    if seqfile.startswith("NC"): #skip salmonella sequences.
   #        gbk_id = seqfile.split(".")[0]
   #        all_nuc_FASTAdicts[gbk_id] = make_FASTAdict('/'.join([full_nuc_path,seqfile]))
   all_prot_FASTAdicts = {}
   full_prot_path = localdir + prot_path
   for seqfile in listdir(full_prot_path):
       if seqfile.startswith("NC"): #skip salmonella sequences.
           gbk_id = seqfile.split(".")[0]
           all_prot_FASTAdicts[gbk_id] = make_FASTAdict('/'.join([full_prot_path,seqfile]))
    
   orthology_file = open(localdir+"data/hatcher-results/pan-9-59.tmp")
   for i, line in enumerate(orthology_file):
       if i == 0:
           continue
       else:
           data = [x.strip() for x in line.split()]
           REL606protein_id =  re.search("NC_012967$lcl|NC_012967\.1_cdsid_(.+?)\t", line).group(1)
           #print REL606protein_id
           ## make nucleotide and protein input FASTA files for this line.
           #nuc_output_dir = localdir + outputdir + "nucleotide/polymorphism/"
           #nuc_output = open(nuc_output_dir + REL606dict[REL606protein_id] + "_ecoli-orthologs.txt", "w")
           prot_output_dir = localdir + outputdir + "protein/polymorphism/"
           prot_output = open(prot_output_dir + REL606dict[REL606protein_id] + "_ecoli-orthologs.txt", "w")
           for elt in data:
               #print elt
               m2 = re.search("(.+?)\$lcl\|(.+?)_cdsid_(.+)", elt)
               gbk = m2.group(1)
               protein_id = m2.group(3)
               #print gbk
               #print protein_id
               #print all_nuc_FASTAdicts[gbk][protein_id]
               #nuc_output.write(all_nuc_FASTAdicts[gbk][protein_id][0])
               #nuc_output.write(all_nuc_FASTAdicts[gbk][protein_id][1])
               #print all_prot_FASTAdicts[gbk][protein_id]
               prot_output.write(all_prot_FASTAdicts[gbk][protein_id][0])
               prot_output.write(all_prot_FASTAdicts[gbk][protein_id][1])
           #nuc_output.close()
           prot_output.close()

def main():
    make_salmonella_ecoli_MAFFT_inputs()
    make_panortholog_MAFFT_inputs()
    

main()
