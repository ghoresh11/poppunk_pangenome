import argparse
import os
import networkx as nx
import csv

## build a massive network where I have
## 39 colors = each of each cluster
## Frequency of each gene in each cluster (size)
## An edge between two genes if in the pairwise roary it was put in the same group
## The COG category of that gene
## Output: Large gene presence absence of all genes
## A network that in an ideal world could be view in cytoscape to see the
## relationships between all the genes

debug = True

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
        if debug and d.split("_")[0] not in ["33","17"]:
                continue
        if os.path.isfile(os.path.join(input_dir, d, "pan_genome_reference.fa")) and os.path.isfile(os.path.join(input_dir, d, "gene_presence_absence.csv")) and os.path.isfile(os.path.join(input_dir, d, "gene_presence_absence.Rtab")):
            dirs_to_return.append(os.path.join(input_dir,d))
    return dirs_to_return


def init_network(orig_roary_dirs):
    ''' initiate all the nodes in the network
    using the original roary output of all the clusters
    Return: the network'''
    print("Initiating properties of all clusters....")
    G = nx.DiGraph()
    cluster_member_gene = {}
    cluster_gene_size = {}
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
                for m in toks[14:]:
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
    ''' add all the COG catgories to the network'''
    print("Getting COG information...")
    cluster_gene_COG = {}

    for cluster in range(1,40):
        if debug and cluster not in [17,33]:
            continue
        cluster = str(cluster)
        cluster_gene_COG[cluster] = {}
        for gene_type in ["rare", "inter", "soft_core", "core"]:
            with open(os.path.join(eggnog_dir, cluster + "_" + gene_type + ".emapper.annotations")) as f:
                for line in f:
                    if line.startswith("#"):
                        continue
                    toks = line.strip().split("\t")
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
            if cluster1 not in ["33","17"] or cluster2 not in ["33","17"]:
                continue

        print("Cluster1: %s, Cluster2: %s" %(cluster1,cluster2))

        gene_members1 = cluster_member_gene[cluster1]
        gene_members2 = cluster_member_gene[cluster2]
        with open(os.path.join(d, "gene_presence_absence.csv")) as f:
            for toks in csv.reader(f):
                flag = False
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
                    else: ## in gene_members2:
                        if gene_members2[m] == "17_cirA_1":
                            flag = True
                        if gene_members2[m] not in cluster2_cnt:
                            cluster2_cnt[gene_members2[m]] = 0
                        cluster2_cnt[gene_members2[m]] += 1

                ### divide the cnts by the size of cluster1 and 2
                # for key in cluster1_cnt:
                #     cluster1_cnt[key] = cluster1_cnt[key] / float(cluster_gene_size[cluster1][key])
                # for key in cluster2_cnt:
                #     cluster2_cnt[key] = cluster2_cnt[key] / float(cluster_gene_size[cluster2][key])

                ## update the edges of the network
                for key in cluster1_cnt:
                    for key2 in cluster2_cnt:
                        G.add_edge(key, key2, proportion=cluster1_cnt[key])
                        G.add_edge(key2, key, proportion=cluster2_cnt[key2])
    return

def generate_cytoscape_output(G, cluster_gene_size, cluster_gene_COG):
    ''' create nodes and edges text files to load in cytoscape'''
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

def merge_genes(orig_roary_dirs, G, cluster_gene_COG):
    ''' use the connected components of the graph to merge
    genes which may not have needed to be seperated in the first place
    and create a new presence-absence file with columns as genes,
    first column is the strain and the second in the cluster (so I can look at
    each cluster indivudially'''
    strains, gene_presence_absence = read_gene_presence_absence(orig_roary_dirs)
    cc = sorted(nx.connected_components(G.to_undirected()), key=len, reverse=True) ## get connected components
    out = open("complete_presence_absence.csv","w")
    out.write("Strain, COG")
    for cluster in range(1,40):
        if debug and cluster not in [17,33]:
            continue
        cluster = str(cluster)
        out.write("," + ",".join(strains[cluster]))
    out.write("\nCluster, COG")
    for cluster in range(1,40):
        if debug and cluster not in [17,33]:
            continue
        cluster = str(cluster)
        out.write("," + ",".join([cluster] * len(strains[cluster])))
    out.write("\n")

    for members in cc:
        cluster_presence_absence = {}
        names = []
        cogs = []
        for m in members:
            cluster = m.split("_")[0]
            names.append("_".join(m.split("_")[1:]))
            if m not in cluster_gene_COG[cluster]:
                cogs.append("?")
            else:
                cogs.append(cluster_gene_COG[cluster][m])
            if cluster not in cluster_presence_absence:
                cluster_presence_absence[cluster] = gene_presence_absence[m]
            cluster_presence_absence[cluster] = [min(sum(x),1) for x in zip(cluster_presence_absence[cluster], gene_presence_absence[m])]
        gene_name = max(set(names), key=names.count)
        cog = max(set(cogs), key=cogs.count)
        out.write(gene_name + "," + cog.replace(",","/"))
        for cluster in range(1,40):
            if debug and cluster not in [17,33]:
                continue
            cluster = str(cluster)
            if cluster not in cluster_presence_absence:
                line = ["0"] * len(strains[cluster])
            else:
                line = map(str,cluster_presence_absence[cluster])
            out.write("," + ",".join(line))
        out.write("\n")
    out.close()
    return cc

def run(args):
    orig_roary_dirs = get_input_dirs(args.input_dir)
    pairwise_roary_dirs = get_input_dirs(args.input_dir_2)
    G, cluster_member_gene, cluster_gene_size = init_network(orig_roary_dirs)
    connect_two_clusters(pairwise_roary_dirs, cluster_member_gene, cluster_gene_size, G)
    cluster_gene_COG = get_COG_cats(cluster_member_gene)
    generate_cytoscape_output(G, cluster_gene_size, cluster_gene_COG)
    cc = merge_genes(orig_roary_dirs, G, cluster_gene_COG)
    analyse_network(G, cc)
    return


def get_options():
    parser = argparse.ArgumentParser(description='Extract the gene sequences from roary outputs, and merge rare genes')
    # input options
    parser.add_argument('--input_dir', required=False,
                        type=str, default =  "/lustre/scratch118/infgen/team216/gh11/e_coli_collections/poppunk/dists_analysis/roary_outputs/",
                        help='path to original roary input directory [%(default)s]')
    parser.add_argument('--input_dir_2', required=False,
                        type=str, default =  "/lustre/scratch118/infgen/team216/gh11/e_coli_collections/poppunk/pairwise_roary/",
                        help='path to pairwise roary directory [%(default)s]')
    return parser.parse_args()


if __name__ == "__main__":
    # get arguments from user
    options = get_options()
    run(options)
