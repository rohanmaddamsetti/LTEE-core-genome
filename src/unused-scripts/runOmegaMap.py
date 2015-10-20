#!/usr/bin/python
##runOmegaMap.py by Rohan Maddamsetti
##This short script is a wrapper for omegaMap.
##A second script will grep through the result files,
##printing out a text file with two headings: locus_tag, and theta_s.

from subprocess import Popen
from os import listdir
import os.path
#from os.path import splitext, basename
import sys
import re

def getLocusTag(path):
    '''This short helper function pulls out a locus_tag from a path.
    '''
    filename = os.path.basename(path)
    pattern = '^(ECB_[0-9]{5})'
#    rest, extension = splitext(path)
    ##pull out the locus_tag from the path.
    try:
        locus_tag = re.search(pattern, filename).group(1)
        #print locus_tag
    except ValueError:
        locus_tag = None
    return locus_tag
    
def createConfigFile (fasta_path, config_path):
    '''This function takes the path to an alignment file,
    and writes a *.ini config file to the second parameter, a directory,
named 'configs'. The * wildcard is replaced by the name of the
fasta file, *.fasta. 

parameter settings for omegaMap (see Luscombe supplement section 2.1.7):
 norders = 10
 niter = 150000
 thinning = 150
 muPrior = improper_inverse
 kappaPrior = improper_inverse
 indelPrior = improper_inverse
 omegaPrior = inverse
 omegaParam = 0.01, 100
 rhoPrior = inverse
 rhoParam = 0.01, 100
 muStart = 0.1
 kappaStart = 3.0
 indelStart = 0.1
 omega_model = variable
 oBlock = 15
 rho_model = variable
 rBlock = 15
    '''

    locus_tag = getLocusTag(fasta_path)

    ##open the init file for writing.
    config_outfile = open(os.path.join(config_path,locus_tag + '.ini'), 'w')

    config_outfile.write("FASTA = " + fasta_path + "\n")
    config_outfile.write("pi = .016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393,.016393" + "\n")
    config_outfile.write("norders = 10" + "\n")
    config_outfile.write("niter = 150000" + "\n")
    config_outfile.write("thinning = 150" + "\n")
    config_outfile.write("muPrior = improper_inverse" + "\n")
    config_outfile.write("kappaPrior = improper_inverse" + "\n")
    config_outfile.write("indelPrior = improper_inverse" + "\n")
    config_outfile.write("omegaPrior = inverse" + "\n")
    config_outfile.write("omegaParam = 0.01, 100" + "\n")
    config_outfile.write("rhoPrior = inverse" + "\n")
    config_outfile.write("rhoParam = 0.01, 100" + "\n")
    config_outfile.write("muStart = 0.1" + "\n")
    config_outfile.write("kappaStart = 3.0" + "\n")
    config_outfile.write("indelStart = 0.1" + "\n")
    config_outfile.write("omega_model = variable" + "\n")
    config_outfile.write("oBlock = 15" + "\n")
    config_outfile.write("rho_model = variable" + "\n")
    config_outfile.write("rBlock = 15" + "\n")

    config_outfile.close()

def createAllConfigs(datapath, config_path):
    '''This function creates all config files for all input files to
    OmegaMap.'''
    datafiles = [x for x in listdir(datapath) if x.startswith("ECB")]
    for base in datafiles:
        path = os.path.join(datapath,base)
        createConfigFile(path, config_path)
    
def runOmegaMap(init_path, resultspath, omegamap_path):
    '''This function runs OmegaMap with the given init file.'''
    name = getLocusTag(init_path)
    out_filename = os.path.join(resultspath, name + ".txt")
    process_args = [omegamap_path, init_path, '-outfile', out_filename]
    #print process_args
    proc =  Popen(process_args, shell=False,
                  stdin=None, stdout=None, stderr=None, close_fds=True)
    return None

def main(argv=None):
    local_path = "/Users/Rohandinho/Desktop/Projects/LTEE-evolutionary-rates/"
    data_path= os.path.join(local_path, "data/omegaMap-input/divergence")
    config_path = os.path.join(local_path, "configs/divergence/")

    createAllConfigs(data_path, config_path)  ##comment out when has already been run.

    #results_path = os.path.join(local_path, "results/omegaMap-output/divergence/")
    #omegaMap_path = os.path.join(local_path, "src/omegamap/omegaMap.dat")
    #test_path = os.path.join(config_path,"ECB_00002.ini")
    #runOmegaMap(test_path, results_path, omegaMap_path)
    #configs = [x for x in listdir("./configs") if getLocusTag(x)]
    #print configs
    
if __name__ == "__main__":
    main()
