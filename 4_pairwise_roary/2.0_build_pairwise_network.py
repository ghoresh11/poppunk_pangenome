import argparse
import os
import networkx as nx
import csv
from Bio.SeqIO.FastaIO import SimpleFastaParser
import subprocess


prefix = "new2_"


def get_input_dirs(cluster1, cluster2, input_dir, pairwise=False):
    ''' check all directories in the input dir
    return a list of directories that have all the required files'''
    print("Getting the input directories...")
    input_dir = os.path.abspath(input_dir)
    directories = os.listdir(input_dir)
    dirs_to_return = []
    curr_prefix = prefix
    if pairwise:
        curr_prefix = ""
    for d in directories:
        flag = False
        if not os.path.isdir(os.path.join(input_dir, d)):
            continue
        if not pairwise and "pairwise" in d:
            continue
        toks = d.split("/")[-1].split("_")
        if len(toks) < 2:
            continue
        first_cluster = toks[0]
        second_cluster = toks[1]
        if pairwise:
            if (first_cluster == cluster1 and second_cluster == cluster2) or (second_cluster == cluster1 and first_cluster == cluster2):
                return(os.path.join(input_dir,d))
        if first_cluster == cluster1 or first_cluster == cluster2:
            dirs_to_return.append(os.path.join(input_dir,d))
            continue
        if second_cluster == cluster1 or second_cluster == cluster2:
            dirs_to_return.append(os.path.join(input_dir,d))
            continue
    return dirs_to_return


def init_network(orig_roary_dirs):
    ''' initiate all the nodes in the network
    using the original roary output of all the two clusters
    Return: the network '''
    print("Initiating properties of all clusters....")
    G = nx.Graph()
    cluster_member_gene = {} ## member of cluster to cluster name
    cluster_gene_size = {} ## cluster name to the size of that cluster
    cluster_gene_sequence = {}
    for d in orig_roary_dirs:
        cluster = d.split("/")[-1].split("_")[0]
        print(cluster)
        cluster_member_gene[cluster] = {}
        cluster_gene_size[cluster] = {}
        cluster_gene_sequence[cluster] = {}
        ## get the members in each gene cluster
        with open(os.path.join(d, prefix + "clustered_proteins")) as f:
            for line in f:
                toks = line.strip().split()
                toks[0] = toks[0].replace(":","")
                num_members = len(toks[1:])
                for m in toks[1:]: ## go over all members of this cluster
                    m = m.split("|")[1]
                    cluster_member_gene[cluster][m] = cluster + "_" + toks[0]
                cluster_gene_size[cluster][cluster + "_" + toks[0]] = num_members
        with open(os.path.join(d, prefix + "pan_genome_reference.fa")) as handle:
            for values in SimpleFastaParser(handle):
                    gene_name = cluster + "_" + values[0].split()[-1]
                    cluster_gene_sequence[cluster][gene_name] = values[1]
        ## get the frequency of each gene cluster and init the graph
        with open(os.path.join(d, prefix + "gene_presence_absence.Rtab")) as f:
            for toks in csv.reader(f, delimiter = "\t"):
                if toks[0] == "Gene":
                    continue
                name = cluster + "_" + toks[0]
                freq = sum(map(int,toks[1:])) / float(len(toks)-1)
                G.add_node(name, freq = freq, cluster=cluster)
    return G, cluster_member_gene, cluster_gene_size, cluster_gene_sequence


def connect_two_clusters(cluster1, cluster2, length_ratio, member_coverage, pairwise_roary_dir, cluster_member_gene, cluster_gene_size, cluster_gene_sequence, G):
    '''using the roary output of two clusters
    add edges where required between two reps of the network
    G is directed and the weight of the edge is the proprtion of members of
    that group that mapped to members in the other group'''
    print("Connecting clusters...")
    gene_members1 = cluster_member_gene[cluster1]
    gene_members2 = cluster_member_gene[cluster2]

    gene_sequences1 = cluster_gene_sequence[cluster1]
    gene_sequences2 = cluster_gene_sequence[cluster2]

    with open(os.path.join(pairwise_roary_dir, "clustered_proteins")) as f: ## new clustered proteins
        for line in f:
            toks = line.strip().split("\t")
            toks[0] = toks[0].replace(":","") # name of current member in new roary
            members = toks[1:] ## get all members of current cluster

            cluster1_cnt = {} ## count how many times members come from patricular gene group in 1 and 2
            cluster2_cnt = {}
            for m in members:
                if m in gene_members1: ## it's in 1
                    if gene_members1[m] not in cluster1_cnt:
                        cluster1_cnt[gene_members1[m]] = 0 ## the key is a gene that should be in the original network
                    cluster1_cnt[gene_members1[m]] += 1
                elif m in gene_members2: ## in gene_members2:
                    if gene_members2[m] not in cluster2_cnt:
                        cluster2_cnt[gene_members2[m]] = 0
                    cluster2_cnt[gene_members2[m]] += 1
                else: ## something weird if happening - it's in neither original roary outut
                    ## This is likely a case where roary did filter a gene
                    ## in one roary run but not another
                    continue

            ### divide the cnts by the size of cluster1 and 2 -> this makes
            ## sense, there's no reason to draw an edge of only 1 of 30 sequences
            # were a match one time...
            for key in cluster1_cnt:
                cluster1_cnt[key] = cluster1_cnt[key] / float(cluster_gene_size[cluster1][key])
            for key in cluster2_cnt:
                cluster2_cnt[key] = cluster2_cnt[key] / float(cluster_gene_size[cluster2][key])

            ## update the edges of the network
            for key in cluster1_cnt:
                for key2 in cluster2_cnt:
                    key1_sequence = float(len(gene_sequences1[key]))
                    key2_sequence = float(len(gene_sequences2[key2]))
                    ratio = min(key1_sequence/key2_sequence, key2_sequence/key1_sequence)
                    ## not enough length coverage or there are not enough members matching (probably due to length)
                    if ratio < length_ratio or cluster1_cnt[key] < member_coverage or cluster2_cnt[key2] < member_coverage: ## apply a more strict threshold -> it's still grouping many
                        continue
                    ## do a pairwise sequence alignment for these two sequences
                    proc = subprocess.Popen(["needleman_wunsch", "--printscores", gene_sequences1[key], gene_sequences2[key2]], stdout=subprocess.PIPE)
                    alignment_score = int(proc.stdout.read().split("\n")[2].split(":")[-1])  # get the alignment score from seq_align output
                    if alignment_score < 100: ## don't combine references that have low alignment scores- something went wrong
                        continue
                    G.add_edge(key, key2, ratio = ratio, aligment_score = alignment_score)
    return


def output_network(cluster1, cluster2, cluster_gene_sequence, G, debug = False):
    ''' use the connected components of the graph to merge
    genes which may not have needed to be seperated in the first place
    and create a new presence-absence file with columns as genes,
    first column is the strain and the second in the cluster (so I can look at
    each cluster indivudially'''
    cc = sorted(nx.connected_components(G), key=len, reverse=True) ##
    out = open(cluster1 + "_" + cluster2 + "_gene_matches.txt", "w")
    for members in cc: ## each members is one gene with all its members
        out.write("\t".join(list(members)) + "\n")
        if debug:
            test = open(cluster1 + "_" + cluster2 + "_check.fa", "w")
            for m in members:
                curr_cluster = m.split("_")[0]
                test.write(">" + m + "\n" + cluster_gene_sequence[curr_cluster][m] + "\n")
            test.close()
            test = open(cluster1 + "_")
            debug = False
    return


def run(args):
    orig_roary_dirs = get_input_dirs(
        args.cluster1, args.cluster2, args.input_dir)  # will not have cluster 21
    pairwise_roary_dir = get_input_dirs(args.cluster1, args.cluster2, args.input_dir_2, pairwise = True)
    G, cluster_member_gene, cluster_gene_size, cluster_gene_sequence = init_network(orig_roary_dirs)
    connect_two_clusters(args.cluster1, args.cluster2, args.length_ratio, args.member_coverage, pairwise_roary_dir, cluster_member_gene, cluster_gene_size, cluster_gene_sequence, G)
    output_network(args.cluster1, args.cluster2, cluster_gene_sequence, G, debug = True)
    return


def get_options():
    parser = argparse.ArgumentParser(
        description='For each pair of clusters, create an output network showing how they combine')
    ''' input options
    to run: write a short .sh script to copy paste to submit 1000 jobs...!
    '''
    parser.add_argument('--cluster1', required=False,
                        type=str, default="10",
                        help='The first cluster to use [%(default)s]')
    parser.add_argument('--cluster2', required=False,
                        type=str, default="20",
                        help='The second cluster to use [%(default)s]')
    parser.add_argument('--length_ratio', required=False,
                        type=float, default="0.8",
                        help='length ratio required to combine two gene clusters [%(default)s]')
    parser.add_argument('--member_coverage', required=False,
                        type=float, default="0.8",
                        help='percent of members that must be shared between two gene clusters to be combined [%(default)s]')
    parser.add_argument('--input_dir', required=False,
                        type=str, default="/lustre/scratch118/infgen/team216/gh11/e_coli_collections/poppunk/new_roary/",
                        help='path to original roary input directory for each cluster [%(default)s]')
    parser.add_argument('--input_dir_2', required=False,
                        type=str, default="/lustre/scratch118/infgen/team216/gh11/e_coli_collections/poppunk/new_roary/pairwise/",
                        help='path to pairwise roary directory [%(default)s]')
    return parser.parse_args()


if __name__ == "__main__":
    # get arguments from user
    options = get_options()
    run(options)
