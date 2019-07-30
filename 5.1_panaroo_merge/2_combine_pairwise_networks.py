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

def split_cluster(H):
    ''' use dbscan to split connected components that don't make sense
    because of spurious matches
    When looking at the plots it's clear that there is structure in these structures
    and that sometimes the merge step incorrectly marks two genes as the same'''
    if len(H) > 500: ## can't deal with such a large network...
        return list(H.nodes())
    H = nx.convert_node_labels_to_integers(H, label_attribute = "orig_name")
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
        if labels[i] != -1 and labels[i] not in members_to_return:
            members_to_return[labels[i]] = []
        members_to_return[labels[i]].append(H.nodes(data=True)[nodes[i]]["orig_name"])

    nx.set_node_attributes(H, label_per_node, 'dbscan')

    ## fix noise: for now just take the closet cluster, might be better to seprate entirely
    for n in H.nodes(data=True):
        if n[1]["dbscan"] != -1:
            continue
        curr_neighbours = H.neighbors(n[0])
        neighbour_clusters = []
        for k in curr_neighbours:
            neighbour_clusters.append(H.nodes[k]['dbscan'])
        curr_cluster = max(set(neighbour_clusters), key = neighbour_clusters.count)
        members_to_return[curr_cluster].append(n["orig_name"])
    return members_to_return


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
        H = G.subgraph(all_members)
        split_members = split_cluster(H)

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
                m = m.split()[0]
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
    quit()
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
    cnt = 0
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
            cnt += 1
            print("Rewriting complete presence absence for gene: %d" %cnt)
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
