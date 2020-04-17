from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.SeqUtils import GC
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.Seq import translate
import os
from numpy import mean,std, median

complete = []
## it failed due to time so running again with appending and skipping what is done (except for very last one)
# complete = set()
# with open("easy_gene_props.csv") as f:
#     for line in f:
#         toks = line.strip().split(",")
#         complete.add(toks[0])
#complete = list(complete)

out = open("easy_gene_props.csv", "w")
out.write("Gene,Property,Mean,Median,Std\n")

gene_cluster_loc = "/nfs/pathogen004/gh11/clusters231019/231019_clusters"
all_genes = os.listdir(gene_cluster_loc)
for gene_file in all_genes:
    if not gene_file.endswith(".fa"):
        continue
    gene = gene_file.split(".")[0]
    if gene in complete and gene != complete[-1]:
        continue
    gene_props = {}
    for key in ["GC","start_codon","length","aromaticity","gravy","instability"]:
        gene_props[key] = []
    curr_fasta_file = os.path.join(gene_cluster_loc, gene_file)
    with open(curr_fasta_file) as handle:
        for values in SimpleFastaParser(handle):
            nucl = values[1]
            prot = translate(values[1])
            gene_props["GC"].append(GC(nucl))
            gene_props["start_codon"].append(nucl.upper()[:3])
            if prot[-1] == "*":
                prot = prot[:-1]
            if "*" in prot:
                continue
            gene_props["length"].append(len(prot))
            if "X" in prot or "J" in prot or "Z" in prot or "B" in prot:
                continue
            props = ProteinAnalysis(prot)
            gene_props["aromaticity"].append(props.aromaticity())
            gene_props["gravy"].append(props.gravy())
            gene_props["instability"].append(props.instability_index())

    out.write(",".join([gene, "num_members", str(len(gene_props["GC"])), "NA"]) + "\n")
    for prop in gene_props:
        if prop == "start_codon":
            num_atg = str(gene_props[prop].count("ATG") / float(len(gene_props[prop])))
            out.write(",".join([gene, prop, num_atg, "NA","NA"]) + "\n")
            continue
        out.write(",".join([gene, prop, str(mean(gene_props[prop])), str(median(gene_props[prop])), str(std(gene_props[prop]))]) + "\n")



out.close()
