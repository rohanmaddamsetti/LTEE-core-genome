#!/usr/bin/env python

##grep-top-mutations3.py by Rohan Maddamsetti
## print line if gene in the line

## pipe in input from stdin.
## for instance:
##  cat ../data/Tenaillon-data/nature18959-s3.csv | python grep-top-mutations2.py

import subprocess
from os import listdir
from os import path
import sys

def main():
    genes = ['actP'
             'arcA',
             'arcB',
             'argR',
             'atoC',
             'atoS',
             'cpxR',
             'crp',
             'envZ',
             'fabF',
             'ftsI',
             'glmS',
             'glnD',
             'glpR',
             'hflB',
             'hsdM',
             'hslU',
             'iap',
             'iclR',
             'infB',
             'infC',
             'lrp',
             'malE',
             'malT',
             'mrdA',
             'mreB',
             'mreC',
             'mreD',
             'nadR',
             'nagC',
             'nusA',
             'plsB',
             'purL',
             'pykA',
             'pykF',
             'queA',
             'rne',
             'rplF',
             'rpoB',
             'rpsD',
             'sapF',
             'smf',
             'spoT',
             'thrA',
             'topA',
             'trkH',
             'yabB',
             'yadG',
             'ybaL',
             'yeeF',
             'yejF',
             'ygjG',
             'yhdG',
             'yicL',
             'yihP',
             'yijC',
             'yqjK']

    for i,line in enumerate(sys.stdin):
        line = line.strip()
        if i == 0:
            print(line)
        else:
            for g in genes:
                if g in line:
                    print(line)

main()
