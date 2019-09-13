import os
import re
import networkx as nx
from Bio.SeqIO.FastaIO import SimpleFastaParser
import subprocess

interpro_scan_folder = "/Users/gh11/poppunk_pangenome/4_pairwise_roary/classify_genes/interpro_scan_results/"
interpro_scan_files = os.listdir(interpro_scan_folder)

gene_classification = {}
for f in interpro_scan_files:
    with open(os.path.join(interpro_scan_folder,f)) as f_open:
        for line in f_open:
            toks = line.strip().split("\t")
            if line.startswith("name"):
                header = toks[1:]
                continue
            gene_classification[toks[0]] = toks[1:]



def get_detailed_functions(filein, fileout):
    all_gene_domains_sigs = {}
    with open(fileout,"w") as out:
        out.write("cluster\tphylogroup\tnum_genes\tgene\t" + "\t".join(header) + "\n")
        with open(filein) as f:
            for line in f:
                if line.startswith("cluster"):
                    continue
                toks = line.strip().split(",")
                genes = toks[2].split("/")
                for g in genes:
                    if g == "":
                        continue
                    g = re.sub(r"[^a-zA-Z0-9_\*-]+", '', g)
                    out.write("\t".join([toks[0],toks[3],toks[1],g] + gene_classification[g]) + "\n")

                    sigs = []
                    for item in gene_classification[g]:
                        all_sigs = item.split("/")[1:]
                        for s in all_sigs:
                            if s not in ["consensus disorder prediction","ATP binding","DNA binding","transmembrane transporter activity","membrane"]:
                                sigs.append(s)

                    all_gene_domains_sigs[g] = {"sigs":sigs, "cluster":toks[0]}
    return all_gene_domains_sigs


def summarise_terms_as_network(all_gene_terms, out):
    G = nx.Graph()
    all_genes = list(all_gene_terms.keys())
    for i in range(0,len(all_genes)-1):
        for j in range(i+1, len(all_genes)):
            geneA = all_genes[i]
            geneB = all_genes[j]
            G.add_node(geneA, desc = "/".join(all_gene_terms[geneA]["sigs"]), cluster = all_gene_terms[geneA]["cluster"])
            G.add_node(geneB, desc = "/".join(all_gene_terms[geneB]["sigs"]), cluster = all_gene_terms[geneB]["cluster"])
            ## get number of terms in common
            intersection = list(set(all_gene_terms[geneA]["sigs"]) & set(all_gene_terms[geneB]["sigs"]))
            if len(intersection) > 1:
                G.add_edge(geneA, geneB, weight = len(intersection), terms = "/".join(intersection))
    nx.write_gml(G, out)
    return G


def get_seqs():
    seqs = {}
    with open("/Users/gh11/poppunk_pangenome/4_pairwise_roary/classify_genes/sequences/Cluster_specific_prot.fa") as handle:
        for values in SimpleFastaParser(handle):
            seqs[values[0].split()[0]] = values[1]
    with open("/Users/gh11/poppunk_pangenome/4_pairwise_roary/classify_genes/sequences/Missing_in_one_prot.fa") as handle:
        for values in SimpleFastaParser(handle):
            seqs[values[0].split()[0]] = values[1]
    with open("/Users/gh11/poppunk_pangenome/4_pairwise_roary/classify_genes/sequences/Ubiquitous__soft_core__prot.fa") as handle:
        for values in SimpleFastaParser(handle):
            seqs[values[0].split()[0]] = values[1]
    return seqs


def generate_connected_component_output(G, outfile, outdir):
    out = open(outfile, "w")
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    out.write("CC,Common_terms,Title," + ",".join(map(str,range(1,52))) + "\n")
    ccs = sorted(nx.connected_components(G), key=len, reverse=True)
    index = 1
    seqs = get_seqs()
    for component in ccs:
        if len(component) < 3:
            break
        ## get a name for this cluster
        names = []
        for node in component:
            names.append(node[:4])
        name = max(set(names), key=names.count) + "_" + str(index)
        index += 1

        ## get the terms in this cluster
        terms = set()
        for edge in G.edges(component, data = True):
            curr_terms = edge[2]["terms"].split("/")
            for t in curr_terms:
                terms.add(t.replace(",","-"))

        fasta_out  = open(os.path.join(outdir, name + ".fa"), "w")

        ## get members of each PopPUNK cluster
        curr_out = [0] * 51
        for node in component:
            curr_cluster = int(G.nodes[node]["cluster"]) - 1
            curr_out[curr_cluster] += 1
            fasta_out.write(">" + node + " (" + str(curr_cluster) + ")" + "\n" + seqs[node] + "\n")
        fasta_out.close()
        out_msa = open(os.path.join(outdir, name + "_msa.fa"), "w")
        p = subprocess.Popen(["mafft", os.path.join(outdir, name + ".fa")], stdout=out_msa, stderr=subprocess.PIPE)
        p.wait()
        out_msa.close()

        out.write(name + "," +  "/".join(terms) + ",TODO," + ",".join(map(str,curr_out)) + "\n")
    out.close()
    return




# #get_detailed_functions("depleted_genes.csv", "depleted_genes_detailed.csv")
# all_specific_terms = get_detailed_functions("specific_genes.csv", "specific_genes_detailed.csv")
# G = summarise_terms_as_network(all_specific_terms, "specific_genes_graph.gml")
# generate_connected_component_output(G, "specific_functional_clusters.csv", "specific")
# ## next step: 1. Get connected components that are larger than x
# ## 2. See which clusters are present in each CC (and how many of each cluster, maybe one has more copies than the rest)
# ## 3. Can i give a very broad description of this cluster based on the shared terms?
# ## maybe later look at the sequences...
#

all_depleted_terms = get_detailed_functions("depleted_genes.csv", "depleted_genes_detailed.csv")
G = summarise_terms_as_network(all_depleted_terms, "depleted_genes_graph.gml")
generate_connected_component_output(G, "depleted_functional_clusters.csv", "depleted")
