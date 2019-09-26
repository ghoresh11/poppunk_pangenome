import os
from numpy import mean,std
import re
''' for the genes that are longer in a cluster, check the length of all the members
of that cluster and the members of the longer cluster to understand if the rep
was chosen wrong, or if it's truly longer'''


clusters_dir = "/nfs/pathogen004/gh11/clusters190919/"
out = open("gene_lengths.csv", "w")
out.write("Group,Gene,Mean,Sd\n")
with open("specific_gene_types.csv") as f:
    for line in f:
        if line.startswith("Gene"):
            continue
        toks = line.strip().split(",")
        if toks[1] not in ["longer","truncated"]:
            continue
        curr_gene = toks[0]
        partners = toks[2].split("/")
        partners.append(curr_gene)
        for gene in partners:
            gene_lengths = []
            if os.path.isfile(os.path.join(clusters_dir, gene + ".fa")):
                with open(os.path.join(clusters_dir, gene + ".fa")) as f2:
                    for line in f2:
                        if not line.startswith(">"):
                            gene_lengths.append(len(line.strip()))
                out.write(curr_gene + "," + gene + "," + str(mean(gene_lengths)) + "," + str(std(gene_lengths)) + "\n")

out.close()
