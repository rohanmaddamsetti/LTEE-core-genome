#!/usr/bin/python

# flatten_Nam_Database.py by Rohan Maddamsetti

# This utility script formats the 'note' entries in rough_Nam_Database.csv.

# Usage: python flatten_Nam_Database.py > ../data/Nam_Database.csv

nam_to_format = open("../data/rough_Nam_Database.csv", "r")

for i, line in enumerate(nam_to_format):
    line = line.strip()
    if i == 0:
        print lineN
    else:
        abbrev, note, enzyme_class  = line.split(',')
        for word in note.split(' '):
            if word.startswith('(') and word.endswith(')'):
                print ','.join([abbrev, word[1:-1], enzyme_class])



