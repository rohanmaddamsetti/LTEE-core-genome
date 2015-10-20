#!/usr/bin/python

## run_summarize.py by Rohan Maddamsetti.
## This wrapper script runs Daniel Wilson's summarize.dat in the OmegaMap package
## on all the output files that I'm analyzing.

from os import listdir
import os.path
import subprocess

#./summarize 50000 genes.out1.txt > genes.summary1.txt 

def main():
    proj_dir = "/Users/Rohandinho/Desktop/Projects/LTEE-evolutionary-rates/"

    summarize_input_dir = os.path.join(proj_dir, "results/omegaMap-output/polymorphism")
    summarize_output_dir = os.path.join(proj_dir, "results/codon-rate-analysis/omegaMap-summaries")
    summarize_path = os.path.join(proj_dir, "src/omegaMap/summarize.dat")
    #print summarize_input_dir
    #print summarize_output_dir
    #print summarize_path
    
    for f in listdir(summarize_input_dir):
        if f.startswith("ECB"):
            #print f
            outfile = os.path.join(summarize_output_dir, "summary-for-" + f)
            #print outfile
            f_path = os.path.join(summarize_input_dir, f)
            my_cmd = ' '.join([summarize_path, "15000", f_path, ">", outfile])
            print my_cmd
            p = subprocess.Popen(my_cmd, shell=True)
            os.waitpid(p.pid, 0)
            
main()
