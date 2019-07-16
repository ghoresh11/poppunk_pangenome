import os
import subprocess
import networkx as nx
import argparse
import operator
import sys
from Bio.SeqIO.FastaIO import SimpleFastaParser


def read_graphs(graph_a, graph_b, name_a, name_b):
    ''' read in the graphs of both runs and merge them into one large graph
    return: 1. The merged graph (no edges between the clusters)
            2. a dictionary converting node ids to centroid name'''

    print("Reading first graph...")
    graphA = nx.read_gml(graph_a)
    for n in graphA.nodes(data=True):
        n[1]["cluster"] = name_a

    print("Reading the second graph...")
    graphB = nx.read_gml(graph_b)
    for n in graphB.nodes(data=True):
        n[1]["cluster"] = name_b

    print("Merging two graphs...")
    G = nx.disjoint_union(graphA, graphB)
    geneid_to_graphid = {}
    graphid_to_geneid = {}
    for n in G.nodes(data=True):
        geneid_to_graphid[n[1]["cluster"] + "|" + n[1]["centroid"]] = n[0]
        graphid_to_geneid[n[0]] = n[1]["cluster"] + "|" + n[1]["centroid"]
    return G, geneid_to_graphid, graphid_to_geneid

def get_cluster_size(gene_presence_absence_file):
    ''' get the size of each cluster to define the minimum weight of an
    edge to use'''
    with open(gene_presence_absence_file) as f:
        for line in f:
            toks = line.strip().split(",")
            index = toks.index("Avg group size nuc")
            return len(toks[index+1:])
    sys.exit("Error: not supposed to have reached this point! Check gene presence absence files!")
    return


def run_cdhit(pan_ref_a, pan_ref_b, name_a, name_b, t, a, removed, geneid_to_graphid):
    ''' merge the two pan reference genome files and run cdhit
    to get new clusters
    return: path to cdhit output'''
    out_ref_tmp = name_a + "_" + name_b + '_pan_genome_reference'
    cdhit_out_tmp = name_a + "_" + name_b + '_cdhit_out'
    cnt = 0
    with open(out_ref_tmp, 'w') as outfile:
        for fname in [pan_ref_a, pan_ref_b]:
            curr = name_a if fname == pan_ref_a else name_b
            with open(fname) as handle:
                for values in SimpleFastaParser(handle):
                    ident = curr + "|" + values[0]
                    if geneid_to_graphid[ident] in removed:
                        continue
                    seq = values[1]
                    if ";" in seq:
                        seq = seq.split(";")[0]
                    outfile.write(">" + ident + "\n" + seq + "\n")
                    cnt += 1
    # -n 10, 11 for thresholds 0.95 ~ 1.0
    # -n 8,9    for thresholds 0.90 ~ 0.95
    # -n 7      for thresholds 0.88 ~ 0.9
    # -n 6      for thresholds 0.85 ~ 0.88
    # -n 5      for thresholds 0.80 ~ 0.85
    # -n 4      for thresholds 0.75 ~ 0.8
    if t > 0.95:
        n = 10
    elif t > 0.9:
        n = 8
    elif t > 0.88:
        n = 7
    elif t > 0.85:
        n = 6
    elif t > 0.8:
        n = 5
    else:
        n = 4

    subprocess.call(["cd-hit-est", "-i", out_ref_tmp, "-o", cdhit_out_tmp,
                     "-c", str(t), "-T", "4", "-d", "0", "-A", str(a), "-n", str(n)])
    os.remove(out_ref_tmp)
    print("Num sequences in the fasta file: %d" %cnt)
    return cdhit_out_tmp

def init_network(geneid_to_graphid, removed):
    ''' init the match network'''
    match_graph = nx.Graph()
    for gene in geneid_to_graphid:
        if geneid_to_graphid[gene] in removed:
            continue
        match_graph.add_node(
            geneid_to_graphid[gene], cluster=gene.split("|")[0])
    return match_graph


def add_edge_to_network(curr_members, match_graph, geneid_to_graphid):
    ''' add an edge to the merged network'''
    keys = list(curr_members.keys())
    for i in range(0, len(keys) - 1):
        for j in range(i + 1, len(keys)):
            first = keys[i]
            second = keys[j]
            if geneid_to_graphid[first] not in match_graph or geneid_to_graphid[second] not in match_graph:
                continue
            if curr_members[first] != curr_members[second]:
                match_graph.add_edge(
                    geneid_to_graphid[first], geneid_to_graphid[second])
    return match_graph


def build_match_network(geneid_to_graphid, cdhit_out, removed):
    ''' build a network matching centroids from
    clusters A and B based on the cdhit output
    return: the matched network'''
    print("Building the match graph from cdhit...")
    match_graph = init_network(geneid_to_graphid, removed)
    num_members = 0
    with open(cdhit_out + ".clstr") as f:
        curr_members = {}
        for line in f:
            if line.startswith(">"):
                match_graph = add_edge_to_network(
                    curr_members, match_graph, geneid_to_graphid)
                num_members += len(curr_members)
                curr_members = {}
                continue
            member = line.strip().split("...")[0].split(">")[-1]
            cluster = member.split("|")[0]
            curr_members[member] = cluster
    # don't forget the last one
    match_graph = add_edge_to_network(
        curr_members, match_graph, geneid_to_graphid)
    num_members += len(curr_members)
    os.remove(cdhit_out)
    os.remove(cdhit_out + ".clstr")
    return match_graph


def sort_neighbours_by_weight(G, source):
    ''' prefer to start with neighbors that have a higher weight '''
    neighbors = list(G.neighbors(source))
    weights = {}
    for n in neighbors:
        curr_weight = int(G.edges[(source, n)]['weight'])
        if curr_weight > 8880:  # not this kind of neighbour
            continue
        weights[n] = curr_weight
    sorted_neighbours = sorted(weights.items(), key=operator.itemgetter(1))
    sorted_neighbours = [i[0] for i in sorted_neighbours]
    return sorted_neighbours


def rec_connect_neighbours(chosen, friend, match_graph, G, removed, threshold, a, min_weights):
    ''' recursively go over the neighbours and find if they
    have a match to each other, in that case we would trust
    they should be together
    additions: added limit to not add an edge if the original edge is below a specific
    threshold'''
    chosen_neighbours = sort_neighbours_by_weight(G, chosen)
    friend_neighbours = sort_neighbours_by_weight(G, friend)
    #print("Going over neighbours of: %s+%s" % (chosen, friend))

    for i in chosen_neighbours:  # over here - when you choose the neighbour choose the best one of chosen
        # already added as a match
        if G.edges[(chosen, i)]['weight'] < min_weights[G.nodes(data=True)[chosen]['cluster']] or G.edges[(chosen, i)]['weight'] > 8880:
            continue
        for j in friend_neighbours:
            # already added as a ,atch
            if G.edges[(friend, j)]['weight'] < min_weights[G.nodes(data=True)[friend]['cluster']] or G.edges[(friend, j)]['weight'] > 8880:
                continue
            if match_graph.has_edge(i, j):
                G.add_edge(i, j, weight=8888, ident=threshold,
                           coverage=a)  # by neighbours
                match_graph.remove_node(i)
                match_graph.remove_node(j)
                removed.append(i)
                removed.append(j)
                G, match_graph, removed = rec_connect_neighbours(
                    i, j, match_graph, G, removed, threshold, a, min_weights)
    return G, match_graph, removed


def choose_nodes(match_graph):
    ''' check if there are nodes with degree 1'''
    degrees = list(match_graph.degree())
    for n in degrees:
        if n[1] == 1:
            chosen = n[0]
            friend = list(match_graph.neighbors(chosen))[0]
            return chosen, friend
    return None, None


def remove_degree_1(chosen, friend, G, match_graph, removed, threshold, a, min_weights):
    # step 1 find a node with degree 1
    #print("Chosen: %s+%s" % (chosen, friend))
    G.add_edge(chosen, friend, weight=9999,
               ident=threshold, coverage=a)  # by degree
    match_graph.remove_node(friend)
    match_graph.remove_node(chosen)
    removed.append(friend)
    removed.append(chosen)
    G, match_graph, removed = rec_connect_neighbours(
        chosen, friend, match_graph, G, removed, threshold, a, min_weights)
    return G, match_graph, removed


def correct_non_1_degree(G, match_graph, name_a, name_b):
    ''' removed edges that have no neighbour support and see if we get new
    one to one matches'''
    cc = nx.connected_components(match_graph)
    edges_to_remove = {}
    has_removed = False
    for c in cc:
        if len(c) == 1:  # no edges
            continue
        cluster_a = {}
        cluster_b = {}
        for member in c:
            if match_graph.nodes(data=True)[member]['cluster'] == name_a:
                cluster_a[member] = list(G.neighbors(member))
            else:
                cluster_b[member] = list(G.neighbors(member))

        for node_a in cluster_a:
            for node_b in cluster_b:
                # same CC but no edge between them
                if not match_graph.has_edge(node_a, node_b):
                    continue
                score = 0
                for n1 in cluster_a[node_a]:
                    for n2 in cluster_b[node_b]:
                        # ignore neighbour from same  cluster
                        if G.nodes(data=True)[n1]['cluster'] == G.nodes(data=True)[n2]['cluster']:
                            continue
                        if G.has_edge(n1, n2):
                            # they share neighbours which match each other
                            score += G[n1][n2]['ident']

                if score == 0:  # remove edges that have no neighbour support
                    match_graph.remove_edge(node_a, node_b)
                    has_removed = True
    nx.write_gml(match_graph, path = "final_match_graph.gml")
    return has_removed, match_graph


def create_merged_output(G, match_graph, name_a, name_b):
    ''' create an an output with a match
    for cluster A to cluster B (keep multiple matches from the match network)'''
    out = open(name_a + "_" + name_b + "_gene_matches.txt", "w")
    singletons = 0
    one_to_one = 0
    one_to_many = 0
    # deal with no one to one cases and singletons
    cc = nx.connected_components(match_graph)
    written = []
    for c in cc:
        members_genes = []
        for m in list(c):
            if m not in written:
                members_genes.append(G.nodes(data=True)[
                                     m]['cluster'] + "_" + G.nodes[m]['name'])
                written.append(m)
        out.write("\t".join(members_genes) + "\n")
        if len(members_genes) > 2:
            one_to_many += 1
        elif len(members_genes) == 2:
            one_to_one += 1
        else:
            singletons += 1

    # step 2: write the matches from G
    for n1 in G.nodes():
        curr_neighbours = list(G.neighbors(n1))
        if n1 in written:
            continue
        out.write(G.nodes(data=True)[n1]
                  ['cluster'] + "_" + G.nodes[n1]['name'])
        written.append(n1)
        neighbour_found = False
        for n2 in curr_neighbours:
            # they're a different cluster
            if G.nodes(data=True)[n1]['cluster'] != G.nodes(data=True)[n2]['cluster']:
                out.write("\t" + G.nodes(data=True)
                          [n2]['cluster'] + "_" + G.nodes[n2]['name'])
                one_to_one += 1
                neighbour_found = True
                written.append(n2)
        out.write("\n")
        if not neighbour_found:
            singletons += 1
    out.close()
    print("Number of one to one matches: %d, Singletons: %d, One to many: %d, Total num genes: %d" %(one_to_one, singletons, one_to_many, one_to_many + one_to_one + singletons))
    return

def run(args):
    G, geneid_to_graphid, graphid_to_geneid = read_graphs(os.path.join(
        args.a, "final_graph.gml"), os.path.join(args.b, "final_graph.gml"),
        args.name_a, args.name_b)
    removed = []
    size_a = get_cluster_size(os.path.join(args.a, "gene_presence_absence.csv"))
    size_b = get_cluster_size(os.path.join(args.b, "gene_presence_absence.csv"))
    min_weights = {args.name_a: size_a * 0.1, args.name_b: size_b * 0.1}
    for alignment_coverage in [0.99, 0.95, 0.9, 0.85, 0.8, 0.75]:
    #for alignment_coverage in [0.9]:
        for threshold in [0.99, 0.95, 0.9, 0.85, 0.8]:
        #for threshold in [0.95, 0.9]:
            print("Running cdhit with threshold %f..." % threshold)
            cdhit_out = run_cdhit(os.path.join(args.a, "pan_genome_reference.fa"), os.path.join(
                args.b, "pan_genome_reference.fa"), args.name_a, args.name_b, threshold, alignment_coverage, removed, geneid_to_graphid)
            match_graph = build_match_network(
                geneid_to_graphid, cdhit_out, removed)
            has_removed = True
            while has_removed:
                chosen, friend = choose_nodes(match_graph)
                while chosen is not None:
                    G, match_graph, removed = remove_degree_1(
                        chosen, friend, G, match_graph, removed, threshold, alignment_coverage, min_weights)
                    chosen, friend = choose_nodes(match_graph)
                has_removed, match_graph = correct_non_1_degree(G, match_graph, args.name_a, args.name_b)

    ## run final time to make sure to add things that couldn't be fixed (or maybe new ones to ones comeup)
    print("Summarising non one-to-one matches left")
    print("coverage: %f, threshold = %f" %(alignment_coverage, threshold))
    cdhit_out = run_cdhit(os.path.join(args.a, "pan_genome_reference.fa"), os.path.join(
        args.b, "pan_genome_reference.fa"), args.name_a, args.name_b, threshold, alignment_coverage, removed, geneid_to_graphid)
    match_graph = build_match_network(geneid_to_graphid, cdhit_out, removed)
    ## to print the remaining network that has more than degree 1:
    nx.write_gml(match_graph, path= args.name_a + "_" + args.name_b + "_remaining_graph.gml")


    print("Writing merged graph to file...")
    nx.write_gml(G, path=args.name_a + "_" + args.name_b + "_graph_merged.gml")
    create_merged_output(G, match_graph, args.name_a, args.name_b)
    ## possible fix: reconnect to genes if makes sense from graph that they should be connected?
    ## I'm not going to apply the fix because I'm going to repeat this many times
    ## so will have a level of confidence when two genes should be in the same group
    return


def get_options():
    parser = argparse.ArgumentParser(
        description='Extract the gene sequences from roary outputs, and merge rare genes')
    # input options
    parser.add_argument('--a', required=False,
                        type=str, default="47/",
                        help='path to panaroo outputs of first group [%(default)s]')
    parser.add_argument('--b', required=False,
                        type=str, default="40/",
                        help='path to panaroo outputs of second group [%(default)s]')
    parser.add_argument('--name_a', required=False,
                        type=str, default="47",
                        help='name to use for cluster a [%(default)s]')
    parser.add_argument('--name_b', required=False,
                        type=str, default="40",
                        help='name to use for cluster a [%(default)s]')
    return parser.parse_args()


if __name__ == "__main__":
    # get arguments from user
    options = get_options()
    run(options)
