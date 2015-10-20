#!/usr/bin/python

##categorize_mobs.py by Rohan Maddamsetti.
##This script counts the different kinds of IS elements in each of the long-term lines.
##The output goes into a csv file called LTEE_40K_IS_events.csv, which is input to
##an R script called transpositionGrapher.R which makes a stacked bar graph.

import csv

all_mobs_input = open("LTEE_40K_IS_events.txt", "r")
output = open("LTEE_40K_IS_events.csv", "w")
myWriter = csv.writer(output)
myWriter.writerow(["Strain", "Event"])

for line in all_mobs_input:
	if line.startswith('#') and not line.startswith("##"):
		strain,sep,rest = line.partition('_40K_')
		strain = strain[1:] #This removes the leading # sign.
	elif line.startswith("MOB"):
		line_data = line.split('\t')
		IStype = line_data[5]
		myWriter.writerow([strain, IStype])
	elif line.startswith("//"):
		strain = None
		
