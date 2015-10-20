#!/usr/bin/python

# radical-AA-exchanges.py by Rohan Maddamsetti.

## IMPORTANT: This script ignores the araA marker mutation so that it isn't spuriously analyzed as a
## nonsynonymous mutation that evolves in parallel.
## The araA marker is at position 70867, and is a T->C D92G mutation.


from __future__ import division
from os import listdir
import os.path
from pprint import pprint
from Bio.Seq import Seq
from Bio import AlignIO

def count_exp_AA_exchanges(alignment_dir):
    expAAs = []
    for aln in listdir(alignment_dir):
#        print aln
        cur_path = os.path.join(alignment_dir, aln)
        try:
            alignment = AlignIO.read(cur_path, "fasta")
            sites = alignment.get_alignment_length()
            for i in range(sites):
                cur_site = str(alignment[:,i])
                first, second = cur_site[0], cur_site[1]
                if (first != "-") and (second!= "-") and (first != second):
                    expAAs.append([first, second])
        except ValueError: ## if MAFFT failed
            continue
    return expAAs
    
def make_dN_list(diff_directory):
    dN_list = []
    for diff in listdir(diff_directory):
        #print diff
        genome = diff.split(".")[0]
        #print genome
        cur_path = os.path.join(diff_directory, diff)
        cur_handle = open(cur_path)
        nonsynonymous_count = 0
        for line in cur_handle:
            ## check if mutation is araA marker--if so, skip.
            if "gene_name=araA" in line and "aa_position=92" in line and "snp_type=nonsynonymous" in line:
                continue
            if line.startswith("SNP") and "snp_type=nonsynonymous" in line:
                nonsynonymous_count = nonsynonymous_count + 1
                #print line
                fields = line.split("\t")
                gene_name, locus_tags, ref_AA, new_AA = None, None, None, None
                for i in fields:
                    if i.startswith("gene_name"):
                        gene_name = i.split("=")[1]
                    if i.startswith("locus_tag"):
                        locus_tags = i.split("=")[1]
                    if i.startswith("aa_ref_seq"):
                        ref_AA = i.split("=")[1]
                    if i.startswith("aa_new_seq"):
                        new_AA = i.split("=")[1]
                if ref_AA and new_AA and new_AA != '*': ## disallow stops.
                    dN_list.append([gene_name, locus_tags, ref_AA, new_AA])
                    #print gene_name, locus_tags, ref_AA, new_AA
    return dN_list

def partition_dNs(dN_list):
    parallels = []
    singletons = []
    par_genes = [] ## for handling a corner case.
    genes = [elt[0] for elt in dN_list]
    #print genes
    for elt in dN_list:
        cur_gene = genes.pop(0)
        if cur_gene in genes or cur_gene in par_genes:
            parallels.append(elt)
            par_genes.append(cur_gene)
        else:
            singletons.append(elt)
    return (parallels, singletons)

## See Dagan and Graur 2002 for these classifications. 

def is_radical_by_charge(pair):
    positive = ['R', 'H', 'K']
    negative = ['D', 'E']
    uncharged = ['A', 'N', 'C', 'Q', 'G', 'I', 'L', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
    i,j = pair
    if i in positive and j in positive:
        return False
    elif i in negative and j in negative:
        return False
    elif i in uncharged and j in uncharged:
        return False
    else:
        return True

def is_radical_by_volume_and_polarity(pair):
    special = ['C']
    neutral_and_small = ['A', 'G', 'P', 'S', 'T'], 
    polar_and_relatively_small = ['N', 'D', 'Q', 'E']
    polar_and_relatively_large = ['R', 'H', 'K']
    nonpolar_and_relatively_small = ['I', 'L', 'M', 'V']
    nonpolar_and_relatively_large = ['F', 'W', 'Y']

    i,j = pair
    if i in neutral_and_small and j in neutral_and_small:
        return False
    elif i in polar_and_relatively_small and j in polar_and_relatively_small:
        return False
    elif i in polar_and_relatively_large and j in polar_and_relatively_large:
        return False
    elif i in nonpolar_and_relatively_small and j in nonpolar_and_relatively_small:
        return False
    elif i in nonpolar_and_relatively_large and j in nonpolar_and_relatively_large:
        return False
    else:
        return True

def is_radical_by_Grantham_Distance(pair):
    radical_pairs = [ ['R','S'],
                      ['L','S'],
                      ['V','S'],
                      ['I','S'],
                      ['F','S'],
                      ['Y','S'],
                      ['C','S'],
                      ['K','S'],
                      ['M','S'],
                      ['W','S'],
                      ['L','R'],
                      ['P','R'],
                      ['A','R'],
                      ['G','R'],
                      ['C','R'],
                      ['W','R'],
                      ['G','L'],
                      ['C','L'],
                      ['Q','L'],
                      ['N','L'],
                      ['K','L'],
                      ['D','L'],
                      ['E','L'],
                      ['F','P'],
                      ['Y','P'],
                      ['C','P'],
                      ['K','P'],
                      ['D','P'],
                      ['W','P'],
                      ['F','T'],
                      ['C','T'],
                      ['W','T'],
                      ['F','A'],
                      ['Y','A'],
                      ['C','A'],
                      ['N','A'],
                      ['K','A'],
                      ['D','A'],
                      ['E','A'],
                      ['W','A'],
                      ['G','V'],
                      ['C','V'],
                      ['N','V'],
                      ['D','V'],
                      ['E','V'],
                      ['I','G'],
                      ['F','G'],
                      ['Y','G'],
                      ['C','G'],
                      ['K','G'],
                      ['M','G'],
                      ['W','G'],
                      ['C','I'],
                      ['Q','I'],
                      ['N','I'],
                      ['K','I'],
                      ['D','I'],
                      ['E','I'],
                      ['C','F'],
                      ['H','F'],
                      ['Q','F'],
                      ['N','F'],
                      ['K','F'],
                      ['D','F'],
                      ['E','F'],
                      ['C','Y'],
                      ['N','Y'],
                      ['D','Y'],
                      ['E','Y'],
                      ['H','C'],
                      ['Q','C'],
                      ['N','C'],
                      ['K','C'],
                      ['D','C'],
                      ['E','C'],
                      ['M','C'],
                      ['W','C'],
                      ['W','H'],
                      ['M','Q'],
                      ['W','Q'],
                      ['M','N'],
                      ['W','N'],
                      ['D','K'],
                      ['W','K'],
                      ['M','D'],
                      ['W','D'],
                      ['M','E'],
                      ['W','E'] ]

    rev_radical_pairs = [ [j,i] for [i,j] in radical_pairs ]
    all_radical_pairs = radical_pairs + rev_radical_pairs
    if pair in all_radical_pairs:
        return True
    else:
        return False

def count_exchanges(l, is_radical_func):
    radicals =  sum([is_radical_func(pair) for pair in l])
    conservatives = len(l) - radicals
    return (radicals, conservatives)

def printContingencyTables(diff_location, alignments_location):
    dN_list = make_dN_list(diff_location)
    parallels, singletons = partition_dNs(dN_list)
    parallel_exchanges = [ x[-2:] for x in parallels ]
    singleton_exchanges = [ x[-2:] for x in singletons ]

    expAA_exchanges = count_exp_AA_exchanges(alignments_location)
    
    print "Radical by Volume and Polarity"
    print "(Radicals, Conservatives)"
    print "Expected:", count_exchanges(expAA_exchanges, is_radical_by_volume_and_polarity)
    print "Parallel:", count_exchanges(parallel_exchanges, is_radical_by_volume_and_polarity)
    print "Singleton:", count_exchanges(singleton_exchanges, is_radical_by_volume_and_polarity)
    print
    print "Radical by Grantham Distance > 100"
    print "(Radicals, Conservatives)"
    print "Expected", count_exchanges(expAA_exchanges, is_radical_by_Grantham_Distance)
    print "Parallel", count_exchanges(parallel_exchanges, is_radical_by_Grantham_Distance)
    print "Singleton", count_exchanges(singleton_exchanges, is_radical_by_Grantham_Distance)


def main():
    project_path = "/Users/Rohandinho/Desktop/Projects/LTEE-core-genome/"
    non_mutator_diff_location = project_path + "/data/annotated_non-mutator_40K_diffs"
    mutator_diff_location = project_path + "/data/annotated_mutator_40K_diffs"
    alignments_location = project_path + "results/MAFFT-output/protein/divergence"

    print "NON-MUTATORS:"
    printContingencyTables(non_mutator_diff_location, alignments_location)
    print "HYPER-MUTATORS:"
    printContingencyTables(mutator_diff_location, alignments_location)

    ## For non-mutators, results depend on how classifying exchanges as radical or conservative.
    ## Test not significant by charge. 
    ## Tests are significant by Volume and Polarity, and by Grantham Distance.

main()
