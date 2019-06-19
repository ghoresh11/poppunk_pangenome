import sys
import os
from Bio.SeqIO.FastaIO import SimpleFastaParser


directory = sys.argv[1]
gene_name = sys.argv[2]

with open(os.path.join(directory, "new_pan_genome_reference.fa")) as handle:
    for values in SimpleFastaParser(handle):
        curr_gene_name = values[0].split()[-1][:-2]
        if curr_gene_name == gene_name:
            print(">" + values[0] + " " + str(len(values[1])) + "\n" + values[1])
