import argparse
import os
import networkx as nx
import csv
from Bio.SeqIO.FastaIO import SimpleFastaParser
import subprocess


files = os.list(".")

## suffix:

G = nx.Graph()
for f in files:
    if not f.endswith("_gene_matches.txt"):
        continue
    with open(f) as f_open:
        for line in f_open:
            members = line.strip().split()
            for m in members: ## make sure you're adding all the clusters once
                G.add_node(m)
            ## add an edge between two genes in the same row of a file
            for i in range(0,len(members)-1):
                for j in range(i+1, len(members)):
                    G.add_edge(members[i], members[j])


## get connected components

## update the files like I did before and merge eveything together
## create the same files as before + a pan_genome_reference file!!
