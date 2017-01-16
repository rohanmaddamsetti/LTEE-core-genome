#!/usr/env python

## process_shiny_data.py by Rohan Maddamsetti

## some python parsing of the Shiny LTEE data before use in plot-KO-parallelism.R

## Usage: python process_shiny_data.py > ../results/processed_shiny_data.csv

from collections import defaultdict
import csv

def print_lineages_per_gene(fname_list):

    ## Shiny CSVs have 21 fields.
    ## 1 row
    ## 2 treatment
    ## 3 type
    ## 4 start_position
    ## 5 end_position
    ## 6 gene_position
    ## 7 html_position
    ## 8 html_mutation
    ## 9 html_mutation_annotation
    ## 10 gene_list
    ## 11 gene_name
    ## 12 html_gene_name
    ## 13 gene_product
    ## 14 html_gene_product
    ## 15 locus_tag
    ## 16 mutation_category
    ## 17 snp_type
    ## 18 num
    ## 19 population
    ## 20 time
    ## 21 strain

    gene_lineage_dict = defaultdict(lambda:'')
    for fname in fname_list:
        with open(fname, 'r') as f:
            for i,row in enumerate(csv.reader(f, delimiter=',', skipinitialspace=True)):
                if i == 0:
                    continue
                if 'intergenic' in row:
                    continue
                if 'synonymous' in row:
                    continue
                gene_list = row[9]
                html_lineages = row[18]
                for g in gene_list.split(','):
                    if g.startswith('[') and g.endswith(']'):
                        g = g[1:-1]
                    if gene_lineage_dict[g] == '':
                        gene_lineage_dict[g] = html_lineages
                    else: ## add a break at the end if concatening existing string.
                        gene_lineage_dict[g] = gene_lineage_dict[g] + '<br>'+ html_lineages

    print('gene_name, lineages')
    for k,v in sorted(gene_lineage_dict.items()):
        lineages = ':'.join(sorted(list(set(v.split('<br>')))))
        print(k+','+lineages)


def main():
    non_mut_file = "../data/Shiny-LTEE-data/all-shiny-non-mutator-data.csv"
    mut_file = "../data/Shiny-LTEE-data/all-shiny-mutator-data.csv"

    print_lineages_per_gene([non_mut_file,mut_file])

main()
