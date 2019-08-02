import argparse
import os
import networkx as nx
import csv
import numpy as np
from Bio.SeqIO.FastaIO import SimpleFastaParser
import subprocess
import time
import glob
import signal
from csv import reader
import igraph
from scipy.special import comb
from sklearn.cluster import dbscan
from math import sqrt
import operator


## build a massive network where I have
## 39 colors = each of each cluster
## Frequency of each gene in each cluster (size)
## An edge between two genes if in the pairwise roary it was put in the same group
## Output: Large gene presence absence of all genes
## A network that in an ideal world could be view in cytoscape to see the
## relationships between all the genes


#clusters_to_remove = [21, 43, 49, 50] ## when not debugging, add more clusters here if they're weird, there is no 50
#clusters_to_remove = range(1,49) + [50] ## when debugging, will run on 24 and 39
clusters_to_remove = [21, 43, 49]
prefix = ""

def get_input_dirs(input_dir, pairwise = False):
    ''' check all directories in the input dir
    return a list of directories that have all the required files'''
    print("Getting the input directories...")
    input_dir = os.path.abspath(input_dir)
    directories = os.listdir(input_dir)
    dirs_to_return = []
    dirs_as_ints = []
    curr_prefix = prefix
    if pairwise:
        curr_prefix = ""

    for d in directories:
        if not os.path.isdir(os.path.join(input_dir, d)):
            continue
        if not pairwise and "pairwise" in d:
            continue
        #  for debugging, uncomment:
        if d.split("_")[0] in map(str, clusters_to_remove):
            continue
        ## use the post-processed outputs
        if os.path.isfile(os.path.join(input_dir, d, curr_prefix + "pan_genome_reference.fa")) and os.path.isfile(os.path.join(input_dir, d, curr_prefix + "gene_presence_absence.Rtab")):
            dirs_to_return.append(os.path.join(input_dir,d))

            if not pairwise:
                dirs_as_ints.append(int(d.replace("/","")))

            if not os.path.isfile(os.path.join(input_dir, d, curr_prefix + "clustered_proteins")):
                presence_absence_to_clustered_proteins(os.path.join(input_dir, d, curr_prefix + "gene_presence_absence.csv"),
                os.path.join(input_dir, d, curr_prefix + "clustered_proteins"),
                os.path.join(input_dir, d, curr_prefix + "gene_data.csv"))
    if not pairwise:
        for i in range(1,52):
            if i not in dirs_as_ints:
                clusters_to_remove.append(i)
    return dirs_to_return

def get_panarooid_to_annotid(gene_data_file):
    ''' parse the gene data file to convert panaroo ID to the annot ID'''
    panarooid_to_annotid = {}
    with open(gene_data_file) as f:
        for line in f:
            toks = line.strip().split(",")
            panarooid_to_annotid[toks[2]] = toks[3]
    return panarooid_to_annotid


def get_panarooid_to_genename(clustered_proteins_file):
    panarooid_to_genename = {}
    with open(clustered_proteins_file) as f:
        for line in f:
            toks = line.strip().split()
            gene_name = toks[0].replace(":","")
            for m in toks[1:]:
                panarooid_to_genename[m] = gene_name
    return panarooid_to_genename

def presence_absence_to_clustered_proteins(in_file, out_file, gene_data_file):
    ''' create a clustered proteins file for the panaroo outputs'''
    panarooid_to_annotid = get_panarooid_to_annotid(gene_data_file)
    out = open(out_file, "w")
    with open(in_file) as f:
        for toks in reader(f, delimiter=','):
            if toks[0] == "Gene":
                member_index = toks.index('Avg group size nuc') + 1
                continue
            out.write(toks[0] + ":")
            for ms in toks[member_index:]:
                if ms == "":
                    continue
                ms = ms.split()
                for m in ms:
                    if "refound" in m:
                        out.write("\t" + m)
                    else:
                        out.write("\t" + panarooid_to_annotid[m])
            out.write("\n")
    out.close()
    return


def init_network(orig_roary_dirs):
    ''' initiate all the nodes in the network
    using the original roary output of all the clusters
    Return: the network '''
    print("Initiating properties of all clusters....")
    G = nx.Graph()
    cluster_member_gene = {} ## member of cluster to cluster name
    cluster_gene_size = {} ## cluster name to the size of that cluster
    for d in orig_roary_dirs:
        cluster = d.split("/")[-1].split("_")[0]
        print(cluster)
        cluster_member_gene[cluster] = {}
        cluster_gene_size[cluster] = {}
        ## get the members in each gene cluster
        with open(os.path.join(d, prefix + "clustered_proteins")) as f:
            for line in f:
                toks = line.strip().split()
                toks[0] = toks[0].replace(":","")
                num_members = len(toks[1:])
                for m in toks[1:]: ## go over all members of this cluster
                    cluster_member_gene[cluster][m] = cluster + "_" + toks[0]
                cluster_gene_size[cluster][cluster + "_" + toks[0]] = num_members
        ## get the frequency of each gene cluster and init the graph
        with open(os.path.join(d, prefix + "gene_presence_absence.Rtab")) as f:
            for toks in csv.reader(f, delimiter = "\t"):
                if toks[0] == "Gene":
                    continue
                name = cluster + "_" + toks[0]
                freq = sum(map(int,toks[1:])) / float(len(toks)-1)
                G.add_node(name, freq = freq, cluster=cluster)
    return G, cluster_member_gene, cluster_gene_size


def connect_two_clusters(pairwise_roary_dirs, cluster_member_gene, cluster_gene_size, G):
    '''using the roary output of two clusters
    add edges where required between two reps of the network
    G is directed and the weight of the edge is the proprtion of members of
    that group that mapped to members in the other group'''
    print("Connecting clusters...")
    print("Calculating pairwise comparisons...")
    for d in pairwise_roary_dirs:
        clusters = d.split("/")[-1].split("_")
        cluster1 = clusters[0]
        cluster2 = clusters[1]

        if int(cluster1) in clusters_to_remove or int(cluster2) in clusters_to_remove:
            continue

        print("Cluster1: %s, Cluster2: %s" %(cluster1,cluster2))

        gene_members1 = cluster_member_gene[cluster1]
        gene_members2 = cluster_member_gene[cluster2]

        with open(os.path.join(d, "clustered_proteins")) as f: ## new presence absence
            for line in f:
                toks = line.strip().split("\t")
                toks[0] = toks[0].replace(":","")
                members = toks[1:] ## get all members of current cluster
                cluster1_cnt = {} ## count how many times members come from patricular gene group in 1 and 2
                cluster2_cnt = {}
                for m in members:
                    if m in gene_members1: ## it's in 1
                        if gene_members1[m] not in cluster1_cnt:
                            cluster1_cnt[gene_members1[m]] = 0
                        cluster1_cnt[gene_members1[m]] += 1
                    elif m in gene_members2: ## in gene_members2:
                        if gene_members2[m] not in cluster2_cnt:
                            cluster2_cnt[gene_members2[m]] = 0
                        cluster2_cnt[gene_members2[m]] += 1
                    else: ## something weird if happening
                        ## This is likely a case where roary did filter a gene
                        ## in one roary run but not another
                        continue

                ### divide the cnts by the size of cluster1 and 2 -> this makes
                ## sense, there's no reason to draw an edge of only 1 of 30 sequences
                # were a match one time...
                for key in cluster1_cnt:
                    cluster1_cnt[key] = cluster1_cnt[key] / float(cluster_gene_size[cluster1][key])
                    G.add_node(key, cluster = cluster1)
                for key in cluster2_cnt:
                    cluster2_cnt[key] = cluster2_cnt[key] / float(cluster_gene_size[cluster2][key])
                    G.add_node(key, cluster = cluster2)

                ## update the edges of the network
                for key in cluster1_cnt:
                    for key2 in cluster2_cnt:
                        if cluster1_cnt[key] > 0.8 and cluster2_cnt[key2] > 0.8: ## limit adding edges majority being matched
                            G.add_edge(key, key2)
    return


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


def align_clusters(G, orig_roary_dirs):
    '''
    This step is a general sanity check to remove spurious edges from the roary outputs
    It does an MSA of all the clusters and removes edges that are more than 50 snp apart.
    I'm assuming that if two genes should be connected some other gene in the middle
    will connect them with fewer than 50 SNPs, otherwise they should be treated
    as distinct genes.
    Go over the network and see if some edges need to be
    removed because of spurious roary matches,
    This can be achieve by aligning all the sequences that form the new cluster
    and see if there is a clear structure within the cluster'''
    print("Creating outputs of merged networks...")
    if not os.path.exists("snp_trees"):
        os.makedirs("snp_trees")
    G = G.to_undirected()
    cc = sorted(nx.connected_components(G), key=len, reverse=True) ## get connected components
    print("PRE")
    print(len(list(cc)))
    cnt = 0
    for members in cc: ## each members is one gene with all its members
        ## first check how many different clusters of each cluster this has
        cnt += 1
        #print("CC: %d..." %cnt)
        tmp_out = open("tmp_cluster_fasta.fa", "w")
        members_per_cluster = {}
        for m in members:
            cluster = m.split("_")[0]
            if cluster not in members_per_cluster:
                members_per_cluster[cluster] = []
            members_per_cluster[cluster].append("_".join(m.split("_")[1:]))
        #
        merged_same_cluster = False
        for val in members_per_cluster.values():
            if len(val) > 1:
                merged_same_cluster = True
                break
        if not merged_same_cluster: ## this cluster only has members from a single original cluster
            continue
        ## write to a temp fasta file to see how they merge using cd-hit
    #    print("getting sequences...")
    #    start_time = time.time()
        num_sequences = 0
        count_num_not_found = 0
        for curr_cluster in members_per_cluster:
            for d in orig_roary_dirs:
                if d.split("/")[-1].split("_")[0] == str(curr_cluster):
                    ref_dir = d
                    if prefix == "panaroo_":
                        panarooid_to_genename = get_panarooid_to_genename(os.path.join(ref_dir, prefix + "clustered_proteins"))
                        panarooid_to_annotid = get_panarooid_to_annotid(os.path.join(ref_dir, prefix + "gene_data.csv"))
                    break
            with open(os.path.join(ref_dir, prefix + "pan_genome_reference.fa")) as handle:
                for values in SimpleFastaParser(handle):
                    if prefix == "panaroo_":
                        if "refound" in values[0]:
                            name = values[0]
                        else:
                            if panarooid_to_annotid[values[0]] in panarooid_to_genename:
                                name = panarooid_to_genename[panarooid_to_annotid[values[0]]]
                            else:
                                name = panarooid_to_annotid[values[0]]
                                count_num_not_found += 1
                    else:
                        name = values[0].split()[1]
                    if name in members_per_cluster[curr_cluster]:
                        ## THis output proves how choosing a different reference can really alter the results.
                        tmp_out.write(">" + curr_cluster + "_" + name + "\n" + values[1].split(";")[0] + "\n")
                        num_sequences += 1
        tmp_out.close()
    #    print("--- %s seconds ---" % (time.time() - start_time))
        if num_sequences > 1000: ## deal with the VERY large cluster seperately (run cdhit-est first...)
            os.rename("tmp_cluster_fasta.fa", str(cnt) + "_cluster_fasta.fa")  ## it needs a lot of cores and time
            ## run cdhit to reduce the number of sequences
            #subprocess.call(["cd-hit-est", "-i", str(cnt) + "_cluster_fasta.fa","-o", str(cnt) + "_cluster_fasta.fa", "-c","0.9", "-T", "4", "-d", "0","-A", "0.8", "-n", "5"])
            continue

        if num_sequences == 1: ## nothing to do
            continue

        ## get the all the SNPs between all the sequences by running an MSA
        out_msa = open("tmp_msa_file.fa", "w")
        p = subprocess.Popen(["mafft", "tmp_cluster_fasta.fa"], stdout=out_msa, stderr=subprocess.PIPE)
        p.wait()
        out_msa.close()
        msa = {}
        ## build a network of pairwise SNP distances
        with open("tmp_msa_file.fa") as handle:
            for values in SimpleFastaParser(handle):
                msa[values[0]] = values[1]

        ## check if there was a cdhit run (when the number of reps is very big):
        cdhit_clusters = {}
        if os.path.isfile(str(cnt) + "_cluster_fasta.fa.clstr"):
            curr_members = []
            curr_rep = ""
            with open(str(cnt) + "_cluster_fasta.fa.clstr") as f:
                for line in f:
                    if line.startswith(">"):
                        cdhit_clusters[curr_rep] = curr_members
                        curr_members = []
                        curr_rep= ""
                        continue
                    member = line.strip().split("...")[0].split(">")[-1]
                    if "*" in line:
                        curr_rep = member
                    curr_members.append(member)
            ## don't forget the last one
            cdhit_clusters[curr_rep] = curr_members

        snp_distances = nx.Graph()
        keys = list(msa.keys())
        for i in range(0, len(keys)-1):
            snp_distances.add_node(keys[i])
            if keys[i] in cdhit_clusters:
                same_as_i = cdhit_clusters[keys[i]]
            else:
                same_as_i = [keys[i]]
            for k in same_as_i:
                snp_distances.add_node(k)
                if k != keys[i]:
                    snp_distances.add_edge(keys[i], k, weight = 0)
            for j in range(i+1, len(keys)):
                if keys[j] in cdhit_clusters:
                    same_as_j = cdhit_clusters[keys[j]]
                else:
                    same_as_j = [keys[j]]
                for k in same_as_j:
                    snp_distances.add_node(k)
                    if k != keys[j]:
                        snp_distances.add_edge(keys[j], k, weight = 0)
                mismatches = sum(1 for x,y in zip(msa[keys[i]].lower(), msa[keys[j]].lower()) if x != y)
                len_ignore = sum(1 for x,y in zip(msa[keys[i]], msa[keys[j]]) if x == y and x == "-") ## not part of the pairwise alignment
                max_mismatches = 0.2 * (len(msa[keys[i]]) - len_ignore) ## this is arbitrary but you have to choose something

                #print("alignment len: %d, shared length: %d, allowed mismatches %f" %(len(msa[keys[i]]),(len(msa[keys[i]]) - len_ignore), max_mismatches))
                if mismatches > max_mismatches:
                    for val1 in same_as_i:
                        for val2 in same_as_j:
                            if G.has_edge(val1,val2):
                                G.remove_edge(val1,val2)
                else:
                    for val1 in same_as_i:
                        for val2 in same_as_j:
                            G.add_edge(val1,val2) ## the edge isn't there when it should be
                    snp_distances.add_edge(keys[i], keys[j], weight = mismatches)


        if snp_distances.number_of_edges() < 2:
            continue

        with open(os.path.join("snp_trees", str(cnt) + "_snp_tree.csv"), "w") as out:
            out.write("Node1, Node2, Cluster, Cluster, SNPs\n")
            ccs = (nx.Graph(snp_distances.subgraph(c)) for c in nx.connected_components(snp_distances))
            for subgraph in ccs:
                subgraph = nx.minimum_spanning_tree(subgraph)
                ## this will print out the original clusters
                for e in snp_distances.edges(data=True):
                    cluster1 = e[0].split("_")[0]
                    cluster2 = e[1].split("_")[0]
                    out.write(",".join([e[0], e[1], cluster1, cluster2, str(e[2]["weight"])]) + "\n")

    for f in glob.glob("tmp*"):
        os.remove(f)
    return G

def split_cluster(name, H):
    ''' use dbscan to split connected components that don't make sense
    because of spurious matches
    When looking at the plots it's clear that there is structure in these structures
    and that sometimes the merge step incorrectly marks two genes as the same'''
    if len(H) > 500: ## can't deal with such a large network... split it first with betweenness centrality
        return [list(H.nodes())]
        #return split_large_cluster(name, H)
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

def merge_clusters(orig_roary_dirs, G):
    ''' use the connected components of the graph to merge
    genes which may not have needed to be seperated in the first place
    and create a new presence-absence file with columns as genes,
    first column is the strain and the second in the cluster (so I can look at
    each cluster indivudially'''
    cc = sorted(nx.connected_components(G), key=len, reverse=True) ##
    print("POST")
    print(len(list(cc)))
    strains, gene_presence_absence = read_gene_presence_absence(orig_roary_dirs)
    # This section is for writing the new presence absence file:
    out = open("complete_presence_absence.csv","w")
    members_out = open("members.csv", "w")

    if not os.path.exists("weird"):
        os.makedirs("weird")

    members_out.write("Gene,Members\n")
    out.write("Strain")
    for cluster in range(1,52):
        if cluster in clusters_to_remove:
            continue
        cluster = str(cluster)
        out.write("," + ",".join(strains[cluster])) ## write all the strains of the cluster
    out.write("\nCluster")
    for cluster in range(1,52):
        if cluster in clusters_to_remove:
            continue
        cluster = str(cluster)
        out.write("," + ",".join([cluster] * len(strains[cluster]))) ## write which cluster they belong to
    out.write("\n")

    used_names = set()
    ## here I can also apply a dbscan step-> similar to what I did for the panaroo outputs...
    for all_members in cc:  # each members is one gene with all its members
        H = nx.Graph(G.subgraph(all_members))
        split_members = split_cluster(get_name_for_cluster(all_members, used_names), H)
        for members in split_members:
            ## first check how many different clusters of each cluster this has
            cluster_counts = {}
            for m in members:
                cluster = m.split("_")[0]
                if cluster not in cluster_counts:
                    cluster_counts[cluster] = 0
                cluster_counts[cluster] += 1
            cluster_presence_absence = {}

            ## init counter for this gene for all clusters
            for cluster in range(1,52):
                if cluster in clusters_to_remove:
                    continue
                cluster_presence_absence[str(cluster)] = [0] * len(strains[str(cluster)])

            names = []
            for m in members:
                cluster = m.split("_")[0] # get the cluster number of this member
                names.append("_".join(m.split("_")[1:])) ## its name
                m = m.split()[0]
                cluster_presence_absence[cluster] = [min(sum(x),1) for x in zip(cluster_presence_absence[cluster], gene_presence_absence[m])] # merge with existing (nothing if the first time)

            gene_name = max(set(names), key=names.count) # get the most common gene name
            while gene_name in used_names: ## add stars to prevent duplicate names
                gene_name = gene_name + "*"

            used_names.add(gene_name)
            members_out.write(gene_name + "," + "\t".join(members) + "\n")

            out.write(gene_name)
            for cluster in range(1,52):
                if cluster in clusters_to_remove:
                    continue
                line = map(str,cluster_presence_absence[str(cluster)])
                out.write("," + ",".join(line))
            out.write("\n")
    out.close()
    members_out.close()
    return

def generate_R_output():
    ''' reopen the complete presence absence file and write files
    that can easily be plotted in R'''
    print("Generating R outputs...")
    melted_freqs = open("melted_gene_freqs.csv","w")
    freqs = open("freqs.csv","w")
    melted_freqs.write("Gene,Cluster,Freq,Class\n")

    freqs.write("Gene")
    for c in range(1,52):
        if c in clusters_to_remove:
            continue
        freqs.write("," + str(c))
    freqs.write("\n")
    with open("complete_presence_absence.csv") as f:
        for line in f:
            if line.startswith("Strain"):
                continue
            toks = line.strip().split(",")
            if line.startswith("Cluster"):
                indices = {}
                for cluster in range(1,52):
                    if cluster in clusters_to_remove:
                        continue
                    indices[cluster] = [i for i, x in enumerate(toks) if x == str(cluster)]
                continue
            gene_name = toks[0]
            curr_freqs = {}
            freqs.write(gene_name)
            for cluster in range(1,52):
                if cluster in clusters_to_remove:
                    continue
                if len(indices[cluster]) == 0:
                    freq = 0
                else:
                    toks = np.array(toks)
                    freq = sum(map(int, toks[indices[cluster]])) / float(len(indices[cluster]))
                if freq < 0.15:
                    gene_class = "rare"
                elif freq < 0.95:
                    gene_class = "inter"
                elif freq < 0.99:
                    gene_class = "soft_core"
                else:
                    gene_class = "core"
                melted_freqs.write(",".join([gene_name, str(cluster), str(freq), gene_class]) + "\n")
                freqs.write(',' + str(freq))
            freqs.write("\n")
    freqs.close()
    melted_freqs.close()
    return

def run(args):
    orig_roary_dirs = get_input_dirs(args.input_dir) ## will not have cluster 21
    pairwise_roary_dirs = get_input_dirs(args.input_dir_2, pairwise = True)
    G, cluster_member_gene, cluster_gene_size = init_network(orig_roary_dirs)
    connect_two_clusters(pairwise_roary_dirs, cluster_member_gene, cluster_gene_size, G)
    G = align_clusters(G, orig_roary_dirs)
    merge_clusters(orig_roary_dirs, G)
    generate_R_output()
    return

def get_options():
    parser = argparse.ArgumentParser(description='Extract the gene sequences from roary outputs, and merge rare genes')
    ''' input options
    to run:
    conda activate python27
    job_name=combine_roary
    bsub -J ${job_name} -R"select[mem>15000] rusage[mem=15000]" -M15000  -G team216 -o ${job_name}.o -e ${job_name}.e python 2_build_pairwise_connections.py
    '''
    parser.add_argument('--input_dir', required=False,
                        type=str, default =  "/lustre/scratch118/infgen/team216/gh11/e_coli_collections/poppunk/new_roary/",
                        help='path to original roary input directory for each cluster [%(default)s]')
    parser.add_argument('--input_dir_2', required=False,
                        type=str, default =  "/lustre/scratch118/infgen/team216/gh11/e_coli_collections/poppunk/new_roary/pairwise/",
                        help='path to pairwise roary directory [%(default)s]')
    return parser.parse_args()


if __name__ == "__main__":
    # get arguments from user
    options = get_options()
    run(options)
