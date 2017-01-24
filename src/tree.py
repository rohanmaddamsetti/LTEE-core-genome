## tree.py by Rohan Maddamsetti.

## 1) make uncorrected distance matrix of dN in core genes
## b/t each pair of strains. Write distance matrix to file.
## 2) Make Neighbor-Joining tree using the distance matrix,
## and write tree to file.
## 3) count dN in each LTEE pop's GD file, ONLY counting
## mutations in core genes. return a dict of distances from REL606.
## 4) Add LTEE pops as children to REL606 node in the NJ tree.
## 5) print tree, and use tree.R to make a tree figure.

import sys
from os import listdir
from os.path import join, exists
import re
import numpy as np
import pandas as pd
from Bio import AlignIO
import pickle
import operator
import itertools
from copy import deepcopy
from skbio import DistanceMatrix
from skbio.tree import nj
from skbio import TreeNode
from ete3 import Tree, TreeStyle, faces, NodeStyle, AttrFace

# compute hamming distance of two sequences
def hamming(s1,s2):
    return sum(map(operator.ne,s1,s2))

def make_distance_matrix(alndir, gbks):
    distmatrix = pd.DataFrame.from_items([(g,[0 for g in gbks]) for g in gbks],
                                         orient='index', columns=gbks)

    for in_aln in [x for x in listdir(alndir) if x.startswith('ECB')]:
        if in_aln.startswith('ECB'):
            locus_tag = in_aln.partition('.')[0]
            try:
                alignment = AlignIO.read(join(alndir,in_aln), "fasta")
                n = len(alignment)
                for i,j in itertools.combinations(range(n),2):
                        sub_aln = alignment[i:j+1:j-i]
                        top = sub_aln[0]
                        bottom = sub_aln[-1]
                        pname = top.id
                        qname = bottom.id
                        p = pname[4:4+9]
                        q = qname[4:4+9]
                        distmatrix[p][q] += hamming(top.seq,bottom.seq)
            except ValueError: ## if MAFFT failed
                continue
    return distmatrix

def count_gd_treecore_dN(treegenes,gd_handle):
    counts = 0
    for l in gd_handle:
        if l.startswith('SNP') and 'nonsynonymous' in l:
            for g in treegenes:
                if g in l:
                    counts += 1
    return counts

def count_LTEE_distance():
    ## only count dN mutations in core genes used in phylogeny construction.
    localdir = "/Users/Rohandinho/Desktop/Projects/LTEE-core-genome/"
    nonmutator_dir = join(localdir,"data/annotated_non-mutator_50K_diffs")
    mutator_dir = join(localdir,"data/annotated_mutator_50K_diffs")
    tree_core_file = open(join(localdir,"results/tree/tree-orthology-summary.csv"))
    tree_core_genes = []
    for i,l in enumerate(tree_core_file):
        if i == 0:
            continue
        l = l.strip()
        ldata = l.split(',')
        tree_core_genes.append(ldata[0].strip(','))

    LTEE_dN_dict = {}
    for f in [x for x in listdir(join(localdir,nonmutator_dir)) if x.endswith('.gd')]:
        lineage = f.split('_')[1]
        fh = open(join(localdir,nonmutator_dir,f))
        LTEE_dN_dict[lineage] = count_gd_treecore_dN(tree_core_genes,fh)
    for f in [x for x in listdir(join(localdir,mutator_dir)) if x.endswith('.gd')]:
        lineage = f.split('_')[1]
        fh = open(join(localdir,mutator_dir,f))
        LTEE_dN_dict[lineage] = count_gd_treecore_dN(tree_core_genes,fh)
    return LTEE_dN_dict


def main():
    localdir = "/Users/Rohandinho/Desktop/Projects/LTEE-core-genome/"
    aln_dir = "results/MAFFT-output/protein/tree/"

    Gs = [ "NC_003197", ## keep the same order as alignments, so lower triangular.
             "NC_000913", "NC_002655", "NC_002695", "NC_004431", "NC_007779", "NC_007946",
             "NC_008253", "NC_008563", "NC_009800", "NC_009801", "NC_010468", "NC_010473",
             "NC_010498", "NC_011353", "NC_011415", "NC_011601", "NC_011741", "NC_011742",
             "NC_011745", "NC_011748", "NC_011750", "NC_011751", "NC_011993", "NC_012759",
             "NC_012892", "NC_012947", "NC_012967", "NC_012971", "NC_013008", "NC_013353",
             "NC_013361", "NC_013364", "NC_013654", "NC_013941", "NC_016902", "NC_017625",
             "NC_017626", "NC_017628", "NC_017631", "NC_017632", "NC_017633", "NC_017634",
             "NC_017635", "NC_017638", "NC_017641", "NC_017644", "NC_017646", "NC_017651",
             "NC_017652", "NC_017656", "NC_017660", "NC_017663", "NC_017664", "NC_017906",
             "NC_018650", "NC_018658", "NC_018661", "NC_020163", "NC_020518", "NC_022364"]

    if not exists(join(localdir,"results/tree/distance_matrix.p")):
        dmatrix = make_distance_matrix(join(localdir,aln_dir),Gs)
        print(dmatrix)
        dmatrix.to_pickle(join(localdir,"results/tree/distance_matrix.p"))
    else:
        dmatrix = pd.read_pickle(join(localdir,"results/tree/distance_matrix.p"))

    ## do some finagling to use NJ algorithm in scikit-bio.
    ## first turn into a numpy array, and add its transpose.
    dist_array = dmatrix.values
    tr_dist_array = np.transpose(deepcopy(dist_array))
    sym_dist_matrix = np.add(dist_array,tr_dist_array)

    ## get the proper names of each bacterial strain.
    ecoli_collection = pd.read_csv(join(localdir,"results/figures-and-supplement/Supp-Table-1-core-genome-strains.txt"),comment='#').drop('lifestyle', 1)
    salmonella = pd.DataFrame(data={'NCBI Refseq ID':'NC_003197','Strain name':"Salmonella enterica LT2"},index=np.arange(1))
    NCBI_data = pd.concat([salmonella,ecoli_collection])
    bacteria_names = NCBI_data['Strain name'].tolist()
    ## annoying hack to write newick tree.
    bacteria_names = [x.replace(' ','_') for x in bacteria_names]

    sk_dmatrix = DistanceMatrix(sym_dist_matrix,bacteria_names)

    ## make the tree with the Neighbor-Joining algorithm.
    sk_tree = nj(sk_dmatrix)

    ## convert skbio tree to ete3 tree.
    tree = Tree.from_skbio(sk_tree)

    ## count dN in the LTEE, only considering tree core genes.
    LTEE_dN_dict = count_LTEE_distance()
    ## add children to the REL606 node.
    REL606_node = tree.search_nodes(name='Escherichia coli B str REL606')[0]
    for k,v in sorted(LTEE_dN_dict.items()):
        REL606_node.add_child(name=k,dist=v)
    ## root the tree.
    root_node = tree.search_nodes(name='Salmonella enterica LT2')[0]
    tree.set_outgroup(root_node)

    ## make a nice tree.

    def layout(node):
        ## modify text by adding attributes to nodes.
        if node.is_leaf():
            faces.add_face_to_node(AttrFace("name",ftype="Arial",text_prefix=" ",fstyle="italic"),node,column=0)

    ts = TreeStyle()
    ts.show_leaf_name = False
    ts.show_branch_length = False
    ts.show_scale=False
    ts.margin_bottom=10
    ts.margin_top=10
    ts.margin_right=10
    ts.layout_fn = layout

    tree.render("/Users/Rohandinho/Desktop/tree.pdf",tree_style=ts)

    print(LTEE_dN_dict)
    ## print distance between REL606 and K-12 as a comparison. (it's 2938.705293.)
    K12_node = tree.search_nodes(name="Escherichia coli K-12 MG1655")[0]
    fundist = tree.get_distance(REL606_node,K12_node)
    print(fundist)

main()
