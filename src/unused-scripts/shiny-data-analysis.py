#!/usr/env python

## shiny-data-analysis.py by Rohan Maddamsetti.

def main():
    mutator_file = open("../data/Shiny-LTEE-data/filtered-shiny-mutators.csv")
    nonmutator_file = open("../data/Shiny-LTEE-data/filtered-shiny-nonmutators.csv")

    for i,line in enumerate(mutator_file):
        line = line.strip()
        if i == 0:
            print(line)
        else:
            fields = line.split(',')



main()
