import os
import re

''' read all the eggnog outputs and summarise COG categories for each of
the 17 gene types'''

## step one-> read in ALL genes from the classification function with
## the "COG" category being "not found"

def get_cat_group(curr_cat):
    if "varied" in curr_cat or curr_cat == "Intermediate_and_rare":
        curr_group = "Varied"
    elif curr_cat.endswith("Cluster_specific"):
        curr_group = "Cluster specific core"
    elif curr_cat in ["Cluster_specific_rare","Multicluster_rare"]:
        curr_group = "Rare"
    elif curr_cat in ["Cluster_specific_intermediate","Multicluster_intermediate"]:
        curr_group = "Intermediate"
    else:
        curr_group = "Core"
    return curr_group


all_genes = {}
with open("gene_classification.csv") as f:
    for line in f:
        toks = line.strip().split("\t")
        curr_cat = re.sub(r"[^a-zA-Z0-9_\*]+", '_', toks[1])
        if "Ubiq" in curr_cat:
            curr_cat = "Ubiquitous"
        all_genes[toks[0]] = curr_cat

out = open("cog_per_gene.csv","w")
out.write("Gene,Group,Cat,COGs\n")
## step two -> read all the eggnog files
## when a gene is given a COG category, remove it from the dictionary
## I can combine the COG cat's with my work on interpro scan for completeness..

eggnog_files = os.listdir("eggnog/")
for f in eggnog_files:
    curr_cat = f.split("_prot.txt")[0]
    curr_cat = re.sub(r"[^a-zA-Z0-9_\*]+", '_', curr_cat)
    if "Ubiq" in curr_cat:
        curr_cat = "Ubiquitous"
    curr_group = get_cat_group(curr_cat)
    with open(os.path.join("eggnog/",f)) as f_open:
        for line in f_open:
            if line.startswith("#"):
                continue
            toks = line.strip().split("\t")
            if len(toks) >= 12:
                curr_cog = toks[-2].replace(",","/").replace(" ","")
                out.write(toks[0] +  "," + curr_group + "," + curr_cat + "," + curr_cog + "\n")
                del all_genes[toks[0]]

for gene in all_genes:
    out.write(gene + ","+ get_cat_group(all_genes[gene]) + "," + all_genes[gene] + ",not_found\n")

out.close()
