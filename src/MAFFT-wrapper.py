#!/usr/bin/python

'''MAFFT-wrapper.py by Rohan Maddamsetti

Uses mafft-linsi to align orthologous proteins among the 60 sampled E. coli
and align orthologous proteins between REL606 and S. typhimurium LT2.

This will be used as input for calculate-divergence.py to:
1) Count the number of sites that have diverged between REL606 and Salmonella 
   in each ortholog.
2) Count the number of polymorphic sites in each panortholog.

'''

from os import listdir, rename, getcwd, remove
import os.path
from sys import stdout, exit
import subprocess
            
def make_alignments(file_of_fasta_sequences, aln_type="divergence", datatype="DNA"):
    '''Input: a FASTA file with the FASTA sequences of gene sequences.
    Output: one file of the aligned genes, one file that contains the gene tree,
    one file of the output, one file of errors encountered, and one file
    that contains the score number of the gene tree.
    '''

    outputdir = "/Users/Rohandinho/Desktop/Projects/LTEE-core-genome/results/MAFFT-output/"
    if datatype == "DNA":
        outputdir = os.path.join(outputdir,"nucleotide",aln_type)
    elif datatype == "Protein":
        outputdir = os.path.join(outputdir,"protein",aln_type)
    filename, sep, ext = os.path.basename(file_of_fasta_sequences).partition('.')
    outfile = os.path.join(outputdir,filename+".aln")
    my_cmd = ' '.join(['mafft-linsi', file_of_fasta_sequences, '>', outfile])
    #print "command is: ", my_cmd
    p = subprocess.Popen(my_cmd, shell=True)
    os.waitpid(p.pid, 0)

def main():

    localdir = "/Users/Rohandinho/Desktop/Projects/LTEE-core-genome/"

    for f in listdir(localdir+"data/MAFFT-input/protein/divergence/"):
    #    print f
        fullpath = localdir+"data/MAFFT-input/protein/divergence/"+f
        make_alignments(fullpath, aln_type="divergence", datatype="Protein")

    for f in listdir(localdir+"data/MAFFT-input/protein/polymorphism/"):
    #    print f
        fullpath = localdir+"data/MAFFT-input/protein/polymorphism/"+f
        make_alignments(fullpath, aln_type="polymorphism", datatype="Protein")
            
if __name__ == '__main__':
    main()

