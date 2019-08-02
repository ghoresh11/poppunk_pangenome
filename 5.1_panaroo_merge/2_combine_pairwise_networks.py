import argparse
import os
import csv
import subprocess
import networkx as nx
import numpy as np
import igraph
from scipy.special import comb
from sklearn.cluster import dbscan
from Bio.SeqIO.FastaIO import SimpleFastaParser
import time
from math import sqrt
import operator

prefix = "panaroo_"


def get_duplicates():
    ''' read in genomes that are duplicated to ignore them in the final output '''
    duplicates = {}
    with open("identical.tab") as f:
        for line in f:
            toks = line.strip().split("\t")
            if toks[0] == toks[1]:
                continue
            name_0 = toks[0].split("/")[-1].split(".")[0]
            name_1 = toks[1].split("/")[-1].split(".")[0]
            if name_0 in duplicates:
                duplicates[name_1] = duplicates[name_0]
            elif name_1 in duplicates:
                duplicates[name_0] = duplicates[name_1]
            else:
                duplicates[name_1] = name_0
                duplicates[name_0] = name_0
    # for d in duplicates:
    #     if duplicates[d] != d:
    #         print(d)
    cnt = {}
    with open("/Users/gh11/e_colis/FILTERED_MD_FINAL_ALL.tab") as f:
        for line in f:
            toks = line.strip().split("\t")
            if line.startswith("ID"):
                gff_loc = toks.index("Annotation_Location")
                poppunk_loc = toks.index("Poppunk_cluster")
                continue
            name = toks[gff_loc].split("/")[-1].split(".")[0]
            if name in duplicates:
                if toks[poppunk_loc] not in cnt:
                    cnt[toks[poppunk_loc]] = 0
                cnt[toks[poppunk_loc]] += 1
    print(poppunk_loc)
    return


def build_graph(combine_dir):
    print("Building the graph")
    files = os.listdir(combine_dir)
    ## suffix:
    G = nx.Graph()
    clusters = set() ## keep a list of clusters that are in the analysis
    for f in files:
        if not f.endswith("_gene_matches.txt"):
            continue
        print(f)
        toks = f.split("_")
        clusters.add(toks[0])
        clusters.add(toks[1])
        with open(f) as f_open:
            for line in f_open:
                members = line.strip().split()
                for m in members: ## make sure you're adding all the clusters once
                    G.add_node(m)
                ## add an edge between two genes in the same row of a file
                for i in range(0,len(members)-1):
                    for j in range(i+1, len(members)):
                        G.add_edge(members[i], members[j])
    return clusters, G

def get_input_dirs(input_dir, clusters):
    ''' check all directories in the input dir
    return a list of directories that have all the required files'''
    print("Getting the input directories...")
    input_dir = os.path.abspath(input_dir)
    directories = os.listdir(input_dir)
    dirs_to_return = []
    curr_prefix = prefix
    for d in directories:
        if not os.path.isdir(os.path.join(input_dir, d)):
            continue
        if "pairwise" in d:
            continue
        #  for debugging, uncomment:
        if d.split("_")[0] not in clusters:
            continue
        ## use the post-processed outputs
        if os.path.isfile(os.path.join(input_dir, d, curr_prefix + "pan_genome_reference.fa")) and os.path.isfile(os.path.join(input_dir, d, curr_prefix + "gene_presence_absence.Rtab")):
            dirs_to_return.append(os.path.join(input_dir,d))
    return dirs_to_return

def read_gene_presence_absence(orig_roary_dirs):
    ''' get the gene presence absence vector for each gene in each
    cluster
    return: a list of strains that will be in the big gene-presence-absence
    file and also the vectors so I can merge vectors of groups that need
    to be merged'''
    print("Reading all the original gene presence absence matrices")
    gene_presence_absence = {}
    strains = {}
    for d in orig_roary_dirs:
        cluster = d.split("/")[-1].split("_")[0]
        with open(os.path.join(d,  prefix + "gene_presence_absence.Rtab")) as f:
            for line in f:
                toks = line.strip().split("\t")
                if line.startswith("Gene"):
                    strains[cluster] = toks[1:]
                    continue
                gene_presence_absence[cluster + "_" + toks[0]] = map(int,toks[1:])
    return strains, gene_presence_absence


def split_large_cluster(name, G):
    ''' when a cluster is too big to handle, split it to a number of
    smaller clusters first based on betweenness_centrality,
    i.e. by removing nodes with very high betweenness_centrality '''
    ## call split_cluster with each subgraph (that way it will carry on until a cluster is size 500)
    ## add the nodes that were removed to which ever cluster they had more edges with
    print("Splitting large cluster using betweenness centrality for size: %d" %len(G))
    k = int(sqrt(len(G)))
    betweenness_centrality = nx.betweenness_centrality(G, k = k)
    betweenness_centrality = sorted(betweenness_centrality.items(), key=operator.itemgetter(1), reverse = True)
    values = [x[1] for x in betweenness_centrality]
    min_val = np.mean(values) + np.std(values)

    remove = {}
    for item in betweenness_centrality:
        if item[1] > min_val:
            remove[item[0]] = set(G.neighbors(item[0]))
            continue
        break
    G.remove_nodes_from(list(remove.keys())) ## remove the nodes with high betweeneness
    ccs_as_list = list(nx.connected_components(G))
    ccs = nx.connected_components(G)
    cnt = 0
    clusters = []
    ## find where the removed nodes should go
    for node in remove:
        G.add_node(node)
        max_neighbours, chosen = -1, -1
        for i in range(0, len(ccs_as_list)):
            num_neighbours = len(list(set(remove[node] & ccs_as_list[i])))
            print("neighbors of removed: %s, CC as a list: %s, Num neighbors: %d" %(str(remove[node]), str(ccs_as_list[i]), num_neighbours))
            if num_neighbours > max_neighbours:
                chosen = i
                max_neighbours = num_neighbours

        ## chosen edges to readd:
        neighbours_to_add = list(set(remove[node] & ccs_as_list[chosen]))
        for nei in neighbours_to_add:
            G.add_edge(node, nei)


    for c in ccs:
        H = nx.Graph(G.subgraph(c))
        clusters += split_cluster(name + "_" + str(cnt), H)
    return clusters


def split_cluster(name, H):
    ''' use dbscan to split connected components that don't make sense
    because of spurious matches
    When looking at the plots it's clear that there is structure in these structures
    and that sometimes the merge step incorrectly marks two genes as the same'''
    if len(H) > 500: ## can't deal with such a large network... split it first with betweenness centrality
        #return [list(H.nodes())]
        return split_large_cluster(name, H)
    if len(H) < 50: ## I'll assume these aren't the problem for now
        return [list(H.nodes())]
    H = nx.convert_node_labels_to_integers(H, label_attribute = "origname")
    edges = zip(*nx.to_edgelist(H))
    H1 = igraph.Graph(len(H), zip(*edges[:2]))
    X = 1 - np.array(H1.similarity_jaccard(loops=False))
    db = dbscan(X, metric='precomputed', eps=0.7, min_samples=4)

    labels = db[1] ## the cluster labeling
    label_per_node = {}
    members_to_return = {}
    nodes = list(H.nodes())
    for i in range(0, len(labels)):
        label_per_node[nodes[i]] = labels[i]
        if labels[i] == -1:
            continue
        if labels[i] not in members_to_return:
            members_to_return[labels[i]] = []
        members_to_return[labels[i]].append(H.nodes[nodes[i]]["origname"])

    nx.set_node_attributes(H, label_per_node, 'dbscan')

    ## fix noise: for now just take the closet cluster, might be better to seprate entirely
    for n in H.nodes(data=True):
        if n[1]["dbscan"] != -1:
            continue
        curr_neighbours = list(H.neighbors(n[0]))
        if len(curr_neighbours) < 3: ## if this node only has 3 edges consider it noise
            print("decided this is noise %s!" %str(curr_neighbours))
            members_to_return[n[1]["origname"]] = [n[1]["origname"]]
            continue
        neighbour_clusters = []
        for k in curr_neighbours:
            neighbour_clusters.append(H.nodes[k]['dbscan'])
        curr_cluster = max(set(neighbour_clusters), key = neighbour_clusters.count)
        if curr_cluster == -1: ## really is noise, it's own cluster
            members_to_return[n[1]["origname"]] = [n[1]["origname"]]
        else:
            members_to_return[curr_cluster].append(n[1]["origname"])
        n[1]["dbscan"] = curr_cluster
    nx.write_gml(H, path = os.path.join("weird", name + ".gml"))
    return list(members_to_return.values())


def get_name_for_cluster(members, used_names):
    names = []
    for m in members:
        cluster = m.split("_")[0] # get the cluster number of this member
        n = "_".join(m.split("_")[1:])
        names.append(n) ## its name
    gene_name = max(set(names), key=names.count) # get the most common gene name
    while gene_name in used_names: ## add stars to prevent duplicate names
        gene_name = gene_name + "*"
    return gene_name

## get connected components
def merge_clusters(orig_roary_dirs, G, clusters):
    ''' use the connected components of the graph to merge
    genes which may not have needed to be seperated in the first place
    and create a new presence-absence file with columns as genes,
    first column is the strain and the second in the cluster (so I can look at
    each cluster indivudially'''
    print("Getting the connected components and creating new files...")
    cc = sorted(nx.connected_components(G), key=len, reverse=True) ##
    strains, gene_presence_absence = read_gene_presence_absence(orig_roary_dirs)
    # This section is for writing the new presence absence file:
    out = open("complete_presence_absence.csv","w")
    members_out = open("members.csv", "w")
    members_out.write("Gene,Members\n")
    out.write("Strain")

    if not os.path.exists("weird"):
        os.makedirs("weird")

    for cluster in range(1,52):
        if str(cluster) not in clusters:
            continue
        cluster = str(cluster)
        out.write("," + ",".join(strains[cluster])) ## write all the strains of the cluster
    out.write("\nCluster")
    for cluster in range(1,52):
        if str(cluster) not in clusters:
            continue
        cluster = str(cluster)
        out.write("," + ",".join([cluster] * len(strains[cluster]))) ## write which cluster they belong to
    out.write("\n")
    weird, good = 0, 0
    used_names = set()
    cnt = 0
    for all_members in cc:  # each members is one gene with all its members
        H = nx.Graph(G.subgraph(all_members))
        split_members = split_cluster(get_name_for_cluster(all_members, used_names), H)
        for members in split_members:
            cnt+=1
            cluster_presence_absence = {}
            ## init counter for this gene for all cluster
            for cluster in range(1,52):
                if str(cluster) not in clusters:
                    continue
                cluster_presence_absence[str(cluster)] = [0] * len(strains[str(cluster)])
            names = []
            for m in members:
                cluster = m.split("_")[0] # get the cluster number of this member
                n = "_".join(m.split("_")[1:])
                names.append(n) ## its name
                cluster_presence_absence[cluster] = [min(sum(x),1) for x in zip(cluster_presence_absence[cluster], gene_presence_absence[m])] # merge with existing (nothing if the first time)

            gene_name = max(set(names), key=names.count) # get the most common gene name
            while gene_name in used_names: ## add stars to prevent duplicate names
                gene_name = gene_name + "*"

            used_names.add(gene_name)
            members_out.write(gene_name + "," + "\t".join(members) + "\n")

            out.write(gene_name)
            for cluster in range(1,52):
                if str(cluster) not in clusters:
                    continue
                line = map(str,cluster_presence_absence[str(cluster)])
                out.write("," + ",".join(line))
            out.write("\n")

    out.close()
    members_out.close()
    print("Number of weird genes: %d, number of good genes: %d" %(weird, good))
    return

def generate_R_output(clusters):
    ''' reopen the complete presence absence file and write files
    that can easily be plotted in R'''
    print("Generating R outputs...")
    freqs = open("freqs.csv","w")
    freqs.write("Gene")
    for c in range(1,52):
        if str(c) not in clusters:
            continue
        freqs.write("," + str(c))
    freqs.write("\n")
    # cnt = 0
    with open("complete_presence_absence.csv") as f:
        for line in f:
            if line.startswith("Strain"):
                continue
            toks = line.strip().split(",")
            if line.startswith("Cluster"):
                indices = {}
                for cluster in range(1,52):
                    if str(cluster) not in clusters:
                        continue
                    indices[cluster] = [i for i, x in enumerate(toks) if x == str(cluster)]
                continue
        #    cnt += 1
        #    print("Rewriting complete presence absence for gene: %d" %cnt)
            gene_name = toks[0]
            freqs.write(gene_name)
            for cluster in range(1,52):
                if str(cluster) not in clusters:
                    continue
                if len(indices[cluster]) == 0:
                    freq = 0
                else:
                    toks = np.array(toks)
                    freq = sum(map(int, toks[indices[cluster]])) / float(len(indices[cluster]))
                freqs.write(',' + str(freq))
            freqs.write("\n")
    freqs.close()
    return


def run():
    clusters, G = build_graph(".")
    input_dirs = get_input_dirs("/lustre/scratch118/infgen/team216/gh11/e_coli_collections/poppunk/new_roary/", clusters)
    merge_clusters(input_dirs, G, clusters)
    generate_R_output(clusters)
    quit()

if __name__ == "__main__":
    run()
