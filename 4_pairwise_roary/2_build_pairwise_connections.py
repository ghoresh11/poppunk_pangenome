import argparse
import os
import networkx as nx
import csv
import numpy as np
from Bio.SeqIO.FastaIO import SimpleFastaParser
from numpy import mean
import subprocess
import time

## build a massive network where I have
## 39 colors = each of each cluster
## Frequency of each gene in each cluster (size)
## An edge between two genes if in the pairwise roary it was put in the same group
## The COG category of that gene
## Output: Large gene presence absence of all genes
## A network that in an ideal world could be view in cytoscape to see the
## relationships between all the genes

debug = True
remove_21 = True # flag to remove cluster 21 as it appears to be clonal

def get_input_dirs(input_dir):
    ''' check all directories in the input dir
    return a list of directories that have all the required files'''
    print("Getting the input directories...")
    input_dir = os.path.abspath(input_dir)
    directories = os.listdir(input_dir)
    dirs_to_return = []
    for d in directories:
        if not os.path.isdir(os.path.join(input_dir, d)):
            continue
        #  for debugging, uncomment:
        if debug and d.split("_")[0] not in ["8","24"]:
                continue
        ## remove cluster 21
        if remove_21 and d.split("_")[0] == "21":
            continue
        if os.path.isfile(os.path.join(input_dir, d, "pan_genome_reference.fa")) and os.path.isfile(os.path.join(input_dir, d, "gene_presence_absence.csv")) and os.path.isfile(os.path.join(input_dir, d, "gene_presence_absence.Rtab")):
            dirs_to_return.append(os.path.join(input_dir,d))
    return dirs_to_return


def init_network(orig_roary_dirs):
    ''' initiate all the nodes in the network
    using the original roary output of all the clusters
    Return: the network '''
    print("Initiating properties of all clusters....")
    G = nx.DiGraph()
    cluster_member_gene = {} ## member of cluster to cluster name
    cluster_gene_size = {} ## cluster name to the size of that cluster
    for d in orig_roary_dirs:
        cluster = d.split("/")[-1].split("_")[0]
        print(cluster)
        cluster_member_gene[cluster] = {}
        cluster_gene_size[cluster] = {}
        ## get the members in each gene cluster
        with open(os.path.join(d, "gene_presence_absence.csv")) as f:
            for toks in csv.reader(f):
                if toks[0] == "Gene":
                    continue
                num_members = 0
                for m in toks[14:]: ## go over all members of this cluster
                    if m == "":
                        continue
                    members = m.split("\t")
                    for m2 in members:
                        num_members += 1
                        cluster_member_gene[cluster][m2] = cluster + "_" + toks[0]
                    cluster_gene_size[cluster][cluster + "_" + toks[0]] = num_members
        ## get the frequency of each gene cluster and init the graph
        with open(os.path.join(d, "gene_presence_absence.Rtab")) as f:
            for toks in csv.reader(f, delimiter = "\t"):
                if toks[0] == "Gene":
                    continue
                name = cluster + "_" + toks[0]
                freq = sum(map(int,toks[1:])) / float(len(toks)-1)
                G.add_node(name, freq = freq, cluster=cluster)
    return G, cluster_member_gene, cluster_gene_size

def get_COG_cats(cluster_member_gene, eggnog_dir = "/lustre/scratch118/infgen/team216/gh11/e_coli_collections/poppunk/eggnog/"):
    ''' add all the COG catgories to the network
    Potential: Skip this and run the EGGNOG analysis only AFTER this step, when genes have been merged'''
    print("Getting COG information...")
    cluster_gene_COG = {}

    for cluster in range(1,40):
        if debug and cluster not in [8,24]:
            continue
        if remove_21 and cluster == 21: ## remove cluster 21
            continue
        cluster = str(cluster)
        cluster_gene_COG[cluster] = {}
        for gene_type in ["rare", "inter", "soft_core", "core"]:
            with open(os.path.join(eggnog_dir, cluster + "_" + gene_type + ".emapper.annotations")) as f:
                for line in f:
                    if line.startswith("#"):
                        continue
                    toks = line.strip().split("\t")
                    if cluster == "1": ## cluster one is special becuase I made the clustered_proteins fasta file
                        gene_name = toks[0]
                    else:
                        gene_name = cluster_member_gene[cluster][toks[0]]
                    if len(toks) < 12:
                        cluster_gene_COG[cluster][gene_name] = "?"
                    else:
                        cluster_gene_COG[cluster][gene_name] = toks[11]
    return cluster_gene_COG


def connect_two_clusters(pairwise_roary_dirs, cluster_member_gene, cluster_gene_size, G):
    '''using the roary output of two clusters
    add edges where required between two reps of the network
    G is directed and the weight of the edge is the proprtion of members of
    that group that mapped to members in the other group'''

    print("Calculating pairwise comparisons...")
    for d in  pairwise_roary_dirs:
        clusters = d.split("/")[-1].split("_")
        cluster1 = clusters[0]
        cluster2 = clusters[1]

        if debug:
            if cluster1 not in ["8","24"] or cluster2 not in ["8","24"]:
                continue

        if remove_21:
            if cluster1 == "21" or cluster2 == "21": ## ignore cluster 21
                continue

        print("Cluster1: %s, Cluster2: %s" %(cluster1,cluster2))

        gene_members1 = cluster_member_gene[cluster1]
        gene_members2 = cluster_member_gene[cluster2]
        with open(os.path.join(d, "gene_presence_absence.csv")) as f: ## new presence absence
            for toks in csv.reader(f):
                if toks[0] == "Gene":
                    continue
                members = [] ## get all members of current cluster
                for m in toks[14:]:
                    if m == "":
                        continue
                    members += m.split("\t")
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
                        ## This is probably a case where roary did filter a gene
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
                            G.add_edge(key, key2, proportion=cluster1_cnt[key])
                            G.add_edge(key2, key, proportion=cluster2_cnt[key2])
    return

def generate_cytoscape_output(G, cluster_gene_size, cluster_gene_COG):
    ''' create nodes and edges text files to load in cytoscape'''
    print("Generating cytoscape output...")
    nodes = open("nodes.txt","w")
    nodes.write("Name,Cluster,Size,Freq, COG\n")
    for n in G.nodes(data=True):
        cluster = n[1]["cluster"]
        gene = n[0]
        cog = "?"
        if gene in cluster_gene_COG[cluster]:
            cog = cluster_gene_COG[cluster][gene]
        nodes.write("\t".join(map(str, [gene, cluster,
        cluster_gene_size[cluster][gene],  n[1]["freq"],
        cog])) + "\n")
    nodes.close()
    edges = open("edges.txt", "w")
    edges.write("Name, Name, Count\n")
    for e in G.edges(data = True):
        edges.write("\t".join(map(str, [e[0], e[1], e[2]["proportion"]])) + "\n")
    edges.close()
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
        with open(os.path.join(d,  "gene_presence_absence.Rtab")) as f:
            for line in f:
                toks = line.strip().split("\t")
                if line.startswith("Gene"):
                    strains[cluster] = toks[1:]
                    continue
                gene_presence_absence[cluster + "_" + toks[0]] = map(int,toks[1:])
    return strains, gene_presence_absence


def remove_edges(G):
    '''Go over the network and see if some edges need to be
    removed because of spurious roary matches'''
    cc = sorted(nx.connected_components(G.to_undirected()), key=len, reverse=True) ## get connected components
    cnt = 0
    for members in cc: ## each members is one gene with all its members
        ## first check how many different clusters of each cluster this has
        tmp_out = open("tmp_cluster_fasta.fa", "w")
        members_per_cluster = {}
        for m in members:
            cluster = m.split("_")[0]
            if cluster not in members_per_cluster:
                members_per_cluster[cluster] = []
            members_per_cluster[cluster].append("_".join(m.split("_")[1:]))

        merged_same_cluster = False
        for val in members_per_cluster.values():
            if len(val) > 1:
                merged_same_cluster = True
                break
        if not merged_same_cluster:
            continue

        expected_num_clusters = mean(map(len, members_per_cluster.values()))
        print(expected_num_clusters)
        ## write to a temp fasta file to see how they merge using cd-hit
        ref_dir = "/lustre/scratch118/infgen/team216/gh11/e_coli_collections/poppunk//dists_analysis/roary_outputs/pan_genome_references/"
        num_sequences = 0
        for curr_cluster in members_per_cluster:
            with open(os.path.join(ref_dir ,curr_cluster + "_reference.fa")) as handle:
                for values in SimpleFastaParser(handle):
                        #print(values[0])
                    if curr_cluster == "1":
                        name = values[0].split()[0]
                    else:
                        name = values[0].split()[1]
                    if name in members_per_cluster[curr_cluster]:
                        tmp_out.write(">" + curr_cluster + "_" + name + "\n" + values[1] + "\n")
                        num_sequences += 1
        tmp_out.close()
        ## get the all the SNPs between all the sequences
        out_msa = open("tmp_msa_file.fa", "w")
        p = subprocess.Popen(["mafft", "--leavegappyregion", "tmp_cluster_fasta.fa"], stdout=out_msa, stderr=subprocess.PIPE)
        p.wait()
        out_msa.close()

        out_snp_sites = open("tmp_snp_sites.fa", "w")
        subprocess.Popen(["snp_sites", "tmp_msa_file.fa"], stdout=out_snp_sites, stderr = subprocess.PIPE)
        out_snp_sites.close()

        snp_distances = [[0] *num_sequences ] * num_sequences
        print(snp_distances)


        cnt += 1
        if cnt == 1:
            quit()

    return


def merge_clusters(orig_roary_dirs, G, cluster_gene_COG):
    ''' use the connected components of the graph to merge
    genes which may not have needed to be seperated in the first place
    and create a new presence-absence file with columns as genes,
    first column is the strain and the second in the cluster (so I can look at
    each cluster indivudially'''

    strains, gene_presence_absence = read_gene_presence_absence(orig_roary_dirs)
    cc = sorted(nx.connected_components(G.to_undirected()), key=len, reverse=True) ## get connected components
    merged_clusters = {}
    ## This section is for writing the new presence absence file:
    # out = open("complete_presence_absence.csv","w")
    #    out.write("Strain, COG")
    # for cluster in range(1,40):
    #     if debug and cluster not in [8,24]:
    #         continue
    #     if remove_21 and cluster == 21: ## remove cluster 21
    #         continue
    #     cluster = str(cluster)
        # out.write("," + ",".join(strains[cluster])) ## write all the strains of the cluster
    # out.write("\nCluster, COG")
    # for cluster in range(1,40):
    #     if debug and cluster not in [8,24]:
    #         continue
    #     if remove_21 and cluster == 21:
    #         continue
    #     cluster = str(cluster)
    #     #out.write("," + ",".join([cluster] * len(strains[cluster])))
    # # out.write("\n")

    used_names = set()
    for members in cc: ## each members is one gene with all its members
        ## first check how many different clusters of each cluster this has
        cluster_counts = {}
        for m in members:
            cluster = m.split("_")[0]
            if cluster not in cluster_counts:
                cluster_counts[cluster] = 0
            cluster_counts[cluster] += 1

        ## if everything is lower than 1, that's fine
        ## if there's even one cluster with 2, we need to check the mean
        ## and possibly separate using cd-hit


        cluster_presence_absence = {}

        ## init counter for this gene for all clusters
        for cluster in range(1,40):
            if debug and cluster not in [8,24]:
                continue
            if remove_21 and cluster == 21:
                continue
            cluster_presence_absence[str(cluster)] = [0] * len(strains[str(cluster)])

        names = []
        cogs = []
        for m in members:
            cluster = m.split("_")[0] # get the cluster number of this member
            names.append("_".join(m.split("_")[1:])) ## its name
        #     if m not in cluster_gene_COG[cluster]:
        #         cogs.append("?")
        #     else:
        #         cogs.append(cluster_gene_COG[cluster][m])
        #     cluster_presence_absence[cluster] = [min(sum(x),1) for x in zip(cluster_presence_absence[cluster], gene_presence_absence[m])] # merge with existing (nothing if the first time)

        gene_name = max(set(names), key=names.count) # get the most common gene name
        while gene_name in used_names: ## add stars to prevent duplicate names
            gene_name = gene_name + "*"

        used_names.add(gene_name)
        merged_clusters[gene_name] = members


        # cog = max(set(cogs), key=cogs.count) ## the COG category is the most common one for this gene
        # out.write(gene_name + "," + cog.replace(",","/"))
        # for cluster in range(1,40):
        #     if debug and cluster not in [8,24]:
        #         continue
        #     if remove_21 and cluster == 21:
        #         continue
        #     line = map(str,cluster_presence_absence[str(cluster)])
        #     out.write("," + ",".join(line))
        # out.write("\n")
    # out.close()
    return merged_clusters

def postprocess_merge(G, merged_clusters):
    ''' Go'''
    for d in orig_roary_dirs:
        clusters.append(d.split("_")[0])
    out = open(argv[2], 'w')
    out.write("Gene, " + ",".join(map(str,clusters)) + ", mean" + "\n")

    with open("tmp_members.csv") as f:
        for line in f:
            if line.startswith("Name"):
                continue
            toks = line.strip().split(",")
            name = toks[0]
            members = toks[1].split("\t")
            counts = {}
            flag = False
            for m in members:
                curr_cluster = int(m.split("_")[0])
                if curr_cluster not in counts:
                    counts[curr_cluster] = 1
                else: ## this cluster was seen already for this gene
                    counts[curr_cluster] += 1
                    flag = True

            if not flag:
                continue

            out.write(name)
            for cluster in clusters:
                if cluster not in counts:
                    out.write(",0")
                else:
                    out.write("," + str(counts[cluster]) )
            out.write("," + str(mean(counts.values())) + "\n")
    out.close()

def generate_R_output():
    ''' reopen the complete presence absence file and write files
    that can easily be plotted in R'''
    print("Generating R outputs...")
    melted_freqs = open("melted_gene_freqs.csv","w")
    freqs = open("freqs.csv","w")
    melted_freqs.write("Gene,Cluster,COG,Freq,Class\n")
    freqs.write("Gene, COG" + ",".join(map(str,range(1,40))) + "\n")
    with open("complete_presence_absence.csv") as f:
        for line in f:
            if line.startswith("Strain"):
                continue
            toks = line.strip().split(",")
            if line.startswith("Cluster"):
                indices = {}
                for cluster in range(1,40):
                    if debug and cluster not in [8,24]:
                        continue
                    if remove_21 and cluster == 21:
                        continue
                    indices[cluster] = [i for i, x in enumerate(toks) if x == str(cluster)]
                continue
            gene_name = toks[0]
            cog = toks[1]
            curr_freqs = {}
            freqs.write(gene_name + "," + cog)
            for cluster in range(1,40):
                if debug and cluster not in [8,24]:
                    continue
                if remove_21 and cluster == 21:
                    continue
                if len(indices[cluster]) == 0:
                    freq = 0
                else:
                    toks = np.array(toks)
                    freq = sum(map(int,toks[indices[cluster]])) / float(len(indices[cluster]))
                if freq < 0.15:
                    gene_class = "rare"
                elif freq < 0.95:
                    gene_class = "inter"
                elif freq < 0.99:
                    gene_class = "soft_core"
                else:
                    gene_class = "core"
                melted_freqs.write(",".join([gene_name, str(cluster), cog, str(freq), gene_class]) + "\n")
                freqs.write("," + str(freq))
            freqs.write("\n")
    freqs.close()
    melted_freqs.close()
    return

def run(args):
    orig_roary_dirs = get_input_dirs(args.input_dir) ## will not have cluster 21
    pairwise_roary_dirs = get_input_dirs(args.input_dir_2)
    G, cluster_member_gene, cluster_gene_size = init_network(orig_roary_dirs)
    cluster_gene_COG = get_COG_cats(cluster_member_gene)
    connect_two_clusters(pairwise_roary_dirs, cluster_member_gene, cluster_gene_size, G)
    generate_cytoscape_output(G, cluster_gene_size, cluster_gene_COG)
    remove_edges(G)
    #merged_clusters = merge_clusters(orig_roary_dirs, G, cluster_gene_COG)
    #postprocess_merge(merged_clusters)
    # generate_R_output()
    return


def get_options():
    parser = argparse.ArgumentParser(description='Extract the gene sequences from roary outputs, and merge rare genes')
    # input options
    ## to run
    # conda activate python27
    # job_name=combine_roary
    # bsub -J ${job_name} -R"select[mem>30000] rusage[mem=30000]" -M30000  -G team216 -o ${job_name}.o -e ${job_name}.e python 2_build_pairwise_connections.py

    parser.add_argument('--input_dir', required=False,
                        type=str, default =  "/lustre/scratch118/infgen/team216/gh11/e_coli_collections/poppunk/dists_analysis/roary_outputs/",
                        help='path to original roary input directory for each cluster [%(default)s]')
    parser.add_argument('--input_dir_2', required=False,
                        type=str, default =  "/lustre/scratch118/infgen/team216/gh11/e_coli_collections/poppunk/pairwise_roary/",
                        help='path to pairwise roary directory [%(default)s]')
    return parser.parse_args()


if __name__ == "__main__":
    # get arguments from user
    options = get_options()
    run(options)
