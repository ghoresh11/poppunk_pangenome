from Bio.SeqIO.FastaIO import SimpleFastaParser

''' get the lengths of ALL the genes which are "truly specific"
Generates a CSV file which can be loaded into R to analyse'''

genes = {}
num_trunc_genes = {}
with open("cluster_all.csv") as f:
    for line in f:
        toks = line.strip().split(",")
        if line.startswith('cluster'):
            cluster_index = toks.index('cluster')
            type_index = toks.index('Type')
            gene_index = toks.index('Gene')
            function_index = toks.index('classification')
            continue

        if toks[cluster_index] not in num_trunc_genes:
            num_trunc_genes[toks[cluster_index]] = 0

        if toks[type_index].strip().lower() == "truncated":
            num_trunc_genes[toks[cluster_index]] += 1
            continue

        if toks[type_index].strip().lower() not in ["true","same functional group"]:
            continue

        functionally_annotated = "No"
        if toks[function_index] != "unknown":
            functionally_annotated = "Yes"
        genes[toks[gene_index]] = {"cluster": toks[cluster_index],
        "functinally_annotated": functionally_annotated }

print("Cluster,Gene,Functionally_annotated,Num_trunc_genes,Length")
with open("../specific.fa") as handle:
    for values in SimpleFastaParser(handle):
        curr_gene = values[0].split()[0]
        if curr_gene not in genes:
            continue
        curr_length = len(values[1])-1
        print(genes[curr_gene]["cluster"] + "," +
            curr_gene + "," +
            genes[curr_gene]["functinally_annotated"] + "," +
            str(num_trunc_genes[genes[curr_gene]["cluster"]]) + "," +
            str(curr_length))
