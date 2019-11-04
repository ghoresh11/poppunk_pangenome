import os
import networkx as nx
import sys

''' go over the interpro-scan results and extract all the genes that I care about
Save those genes, then go over the members file to see which clusters they're in
Generate a network where an edge is the Jaccard similarity between the cluster
each gene has...
I would expect to see that some functions, for instance I can see that is the case with
the autotransporters, they're quite a clear separation between B2 and the rest of the
lineages- but perhaps cluster 40 is more similar than B2??

!!! Outputs !!!
1. keyword_cluster_specific.csv = details all the cluster specific genes and if
they're core/inter/rare
2. keyword_multicluster.csv = can be plotted as a heatmap in the next R script
to look at the distributions of particular functions in the population
3. keyword.gml = a graph that can be loaded into R where there's an edge between
two genes if they have a Jaccard similarity of 0.8 or more. By looking at the
properties of this graph I can see how many times particular functions are
enriched in a lineage (there's a large CC with only one cluster), or there
are many different CCs with only a single cluster (always lineage specific),
or if a function is shared across clusters.
'''

keyword = sys.argv[1].lower()

interproscan_file = "/Users/gh11/poppunk_pangenome/5_classify_genes/interproscan_results.csv"

## get all the genes that have this keyword
all_genes = []
with open(interproscan_file) as f_open:
    for line in f_open:
        if keyword in line.lower():
            all_genes.append(line.strip().split()[0])

## get the phylogroup of all the poppunk clusters
cluster_to_phylo = {}
with open("/Users/gh11/Submissions/my_thesis/Chapter3/figures/cluster_graphics.csv") as f:
    for line in f:
        toks = line.strip().split(",")
        cluster_to_phylo[toks[0]] = toks[1]

## get the classification of all the genes
gene_to_class = {}
with open("/Users/gh11/poppunk_pangenome/5_classify_genes/gene_classification.csv") as f:
    for line in f:
        toks = line.strip().split(",")
        gene_to_class[toks[0]] = toks[1]

## ## create a graph where i keep track of which clusters have this particular gene
nodes = {}
G = nx.Graph()
with open("/Users/gh11/poppunk_pangenome/4_pairwise_roary/231019_corrected/members.csv") as f:
    for line in f:
        toks = line.strip().split(",")
        if toks[0] not in all_genes:
            continue
        nodes[toks[0]] = set()
        phylogroups = set()
        for m in toks[1].split():
            nodes[toks[0]].add(m.split("_")[0])
            phylogroups.add(cluster_to_phylo[m.split("_")[0]])
        G.add_node(toks[0], members = ",".join(list(nodes[toks[0]])), phylogroups = ",".join(list(phylogroups)), geneclass = gene_to_class[toks[0]])


freq_per_node_per_cluster = {}
with open("/Users/gh11/poppunk_pangenome/4_pairwise_roary/231019_corrected/freqs.csv") as f:
        for line in f:
            if line.startswith("Gene"):
                clusters = line.strip().split(",")[1:]
                continue
            toks = line.strip().split(",")
            if toks[0] not in nodes:
                continue
            freqs = toks[1:]
            freq_per_node_per_cluster[toks[0]] = {}
            for i in range(0,len(freqs)):
                freq_per_node_per_cluster[toks[0]][clusters[i]] = float(freqs[i])

def jaccard_similarity(list1, list2):
    list2 = list(list2)
    intersection = list(set(list1).intersection(list2))
    union = (len(list1) + len(list2)) - len(intersection)
    return (float(len(intersection)) / union, ",".join(intersection))


ordered_nodes = list(nodes.keys())
for i in range(0, len(ordered_nodes)-1):
    for j in range(i+1, len(ordered_nodes)):
        score, intersection = jaccard_similarity(nodes[ordered_nodes[i]], nodes[ordered_nodes[j]])
        if score > 0.8:
            G.add_edge(ordered_nodes[i], ordered_nodes[j], score = score, intersection = intersection)


nx.write_gml(G, "functions/" + keyword.replace(" ","_") + ".gml")
## summarise this keyword to plot into R and that way can discuss specificity of
## particular functions
ccs = nx.connected_components(G)
out1 = open("functions/" + keyword.replace(" ","_") + "_cluster_specific.csv", "w")
out1.write("Cluster,Type,Count\n")
out2 = open("functions/" + keyword.replace(" ","_") + "_multicluster.csv", "w")
out2.write("Gene,Cluster,Class,Freq\n")
written_out1 = set()
for c in ccs:
    counts = {"core":0, "rare": 0, "inter": 0}
    clusters = []
    for gene in c:
        clusters += G.node[gene]["members"].split(",")
    clusters = set(clusters)
    if len(clusters) == 1:
        curr_cluster = list(clusters)[0]
        for gene in c:
            if freq_per_node_per_cluster[gene][curr_cluster] < 0.15:
                counts["rare"] += 1
            elif freq_per_node_per_cluster[gene][curr_cluster] < 0.90:
                counts["inter"] += 1
            else:
                counts["core"] += 1
        for t in counts:
            out1.write(curr_cluster + "," + t + "," + str(counts[t]) + "\n")
            written_out1.add(curr_cluster)

    for gene in c:
        for cluster in freq_per_node_per_cluster[gene]:
            out2.write(",".join([gene, cluster, gene_to_class[gene], str(freq_per_node_per_cluster[gene][cluster])]) + "\n")

## add missing clusters to the specific output
for i in range(1,52):
    if i in [21,43,49,50]:
        continue
    if str(i) in written_out1:
        continue
    for t in counts:
        out1.write(str(i) + "," + t + ",0\n")


out1.close()
out2.close()
