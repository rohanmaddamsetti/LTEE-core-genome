#!/usr/bin/python

##simgenome.py by Rohan Maddamsetti.
##This short script makes a probability distribution for where mutations fall
##within protein-coding regions. This is to simulate a null distribution found
##in mutator_analysis.R.

##Go through each CDS, and append the length of the region to a list.
##Sum up the lengths of each CDS examined.
##Then, divide each element of the list by the total sum of nucleotides.
##Output this list to file, for mutator_analysis.R to use as input.

from Bio import SeqIO

cds_lengths = []

for seq_record in SeqIO.parse("../data/REL606.gbk", "genbank"):
    for feature in seq_record.features:
        if feature.type != 'CDS':
            continue
        else:
            cur_start = feature.location.nofuzzy_start
            cur_end = feature.location.nofuzzy_end
            cur_cds_length = cur_end - cur_start
            cds_lengths.append(cur_cds_length)

tot_cds_length = sum(cds_lengths)
normalized_cds_lengths = [float(x)/float(tot_cds_length) for x in cds_lengths]

##output the list of probabilities that a mutation falls in a region to file.
outfile = open("REL606_genebin_probabilities.txt",'w')
probabilities = "\t".join([str(x) for x in normalized_cds_lengths])
outfile.write(probabilities + "\n")
