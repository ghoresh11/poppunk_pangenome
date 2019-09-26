import os
import re
import networkx as nx
from Bio.SeqIO.FastaIO import SimpleFastaParser
import subprocess
import numpy as np
import igraph
from sklearn.cluster import DBSCAN
import operator

def get_gene_classification():
    ''' parse the outputs of interpro scan and ecosyc to get the functional classification
    of all the genes'''
    print("Reading gene classification into dictionary")
    interpro_scan_folder = "/Users/gh11/poppunk_pangenome/5_classify_genes/interpro_scan_results/"
    interpro_scan_files = os.listdir(interpro_scan_folder)
    gene_classification = {}
    for f in interpro_scan_files:
        with open(os.path.join(interpro_scan_folder,f)) as f_open:
            for line in f_open:
                toks = line.strip().split("\t")
                if line.startswith("name"):
                    header = toks[1:]
                    continue
                gene_classification[toks[0]] = ["-"] + toks[1:]
    ## add the ecysyc results
    with open("ecosyc_results.csv") as f:
        for line in f:
            if line.startswith("Gene"):
                continue
            toks = line.strip().split(",")
            gene_classification[toks[0]][0] = toks[1]
    return gene_classification, header

def get_detailed_functions(filein, fileout):
    all_gene_domains_sigs = {}
    with open(fileout,"w") as out:
        out.write("cluster\tphylogroup\tnum_genes\tgene\tpathway\t" + "\t".join(header) + "\n")
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
                    i = 0
                    for item in gene_classification[g]:
                        all_sigs = item.split("/")
                        i += 1
                        if i == 2: ## ignore my classification
                            continue
                        for s in all_sigs:
                            ## ignore some terms that are too general
                            if s not in ["","unknown","DUF","consensus disorder prediction","ATP binding","cytoplasm","catalytic activity",
                            "oxidation-reduction process","integral component of membrane","transporter activity",
                            "DNA binding","transmembrane transporter activity","membrane","-","general function prediction",
                            "Amino acid","outer membrane","P-loop containing nucleoside triphosphate hydrolase","pilus",
                            "ATPase activity"]:
                                sigs.append(s)
                    all_gene_domains_sigs[g] = {"sigs":sigs, "cluster":toks[0]}
    return all_gene_domains_sigs


def summarise_terms_as_network(all_gene_terms, out):
    '''create a network with edges between the nodes if they share
    any functional terms'''
    print("Summarising all the terms in a network...")
    G = nx.Graph()
    all_genes = list(all_gene_terms.keys())
    for i in range(0,len(all_genes)-1):
        for j in range(i, len(all_genes)):
            geneA = all_genes[i]
            geneB = all_genes[j]
            if geneA == geneB: ## ignore self loops
                continue
            G.add_node(geneA, desc = "/".join(all_gene_terms[geneA]["sigs"]), cluster = all_gene_terms[geneA]["cluster"])
            G.add_node(geneB, desc = "/".join(all_gene_terms[geneB]["sigs"]), cluster = all_gene_terms[geneB]["cluster"])
            ## get number of terms in common
            intersection = list(set(all_gene_terms[geneA]["sigs"]) & set(all_gene_terms[geneB]["sigs"]))
            if len(intersection) > 1:
                G.add_edge(geneA, geneB, weight = len(intersection), terms = "/".join(intersection))
    nx.write_gml(G, out)
    return G

def set_specific_type(merged_graph, specific_type, key, new_type, missing = "NA",desc = "-"):
    ''' this is complicated, but there's a hierarchy in defining the property of each
    missing or cluster specific gene, depending on their relationships
    Basically wrong>truncated/missing>same functional group'''
    s_key = "S-" + re.sub(r"[^a-zA-Z0-9_\*-]+", '', key)
    m_key = "M-" + re.sub(r"[^a-zA-Z0-9_\*-]+", '', missing)
    if "property" not in merged_graph.nodes[m_key]:
        merged_graph.nodes[m_key]["property"] = "blast"
        merged_graph.nodes[m_key]["cluster"] = "DB"
        merged_graph.nodes[m_key]["desc"] = "NA"

    orig = specific_type[key]["type"]
    if new_type == "wrong" or specific_type[key]["type"] == "wrong":
        specific_type[key]["type"] = "wrong"
        merged_graph.nodes[s_key]["property"] = "wrong"
        merged_graph.nodes[m_key]["property"] = "wrong"
        specific_type[key]["missing"]["wrong"] = [missing]
        return merged_graph

    if new_type in ["longer","truncated"]:
        specific_type[key]["type"] = new_type
        if new_type not in specific_type[key]["missing"]:
            specific_type[key]["missing"][new_type] = set()
        specific_type[key]["missing"][new_type].add(missing)

    if specific_type[key]["type"] in ["longer","truncated"]:
        merged_graph.nodes[s_key]["property"] = specific_type[key]["type"]
        if merged_graph.nodes[m_key]["property"] != "wrong":
            merged_graph.nodes[m_key]["property"] = specific_type[key]["type"]
        return merged_graph

    if new_type == "same functional group" or specific_type[key]["type"] == "same functional group":
        specific_type[key]["type"] = "same functional group"
        specific_type[key]["desc"] = desc
        merged_graph.nodes[s_key]["property"] =  "same functional group"
        if merged_graph.nodes[m_key]["property"] not in ["longer","truncated","wrong"]:
            merged_graph.nodes[m_key]["property"] =  "same functional group"
        specific_type[key]["missing"]["same functional group"] = [missing]
    return merged_graph


def split_cluster(H, merged_graph):
    ''' use dbscan to split connected components that don't make sense
    because of spurious matches
    When looking at the plots it's clear that there is structure in these structures
    and that sometimes the merge step incorrectly marks two genes as the same'''
    # if len(H) < 20: ## assume these aren't the problem, they only have one member from each cluster?
    #     return [], merged_graph
    H = nx.convert_node_labels_to_integers(H, label_attribute = "origname")
    edges = list(zip(*nx.to_edgelist(H)))
    H1 = igraph.Graph(len(H), list(zip(*edges[:2])))
    X = 1 - np.array(H1.similarity_jaccard(loops=False))
    db = DBSCAN(metric='precomputed', eps=0.5, min_samples=4, n_jobs=16).fit(X)
    labels = db.labels_
    nodes = list(H.nodes())
    for i in range(0, len(labels)):
        H.nodes[nodes[i]]["dbscan"] = labels[i]
        merged_graph.nodes[H.nodes[nodes[i]]["origname"]]["dbscan"] = labels[i]
    ## fix noise
    for n in H.nodes(data=True):
        if n[1]["dbscan"] != -1:
            continue
        curr_neighbours = list(H.neighbors(n[0]))
        if len(curr_neighbours) < 3: ## if this node only has 2 edges consider it noise
            continue
        neighbour_clusters = []
        for k in curr_neighbours:
            neighbour_clusters.append(H.nodes[k]['dbscan'])
        curr_cluster = max(set(neighbour_clusters), key = neighbour_clusters.count)
        if curr_cluster != -1: ## really is noise, it's own cluster
                merged_graph.nodes[n[1]["origname"]]["dbscan"] = curr_cluster
                n[1]["dbscan"] = curr_cluster
    edges_to_remove = []
    for e in H.edges():
        if H.nodes[e[0]]['dbscan'] != H.nodes[e[1]]['dbscan'] or H.nodes[e[1]]['dbscan'] == -1 or H.nodes[e[0]]['dbscan'] == -1:
            edges_to_remove.append((H.nodes[e[0]]['origname'],H.nodes[e[1]]['origname']))
    return edges_to_remove, merged_graph

def connect_both_graphs(specific, depleted):
    ''' in order to understand whether there are genes in a cluster that
    are cluster specific because they are compensated by a different gene with the
    same function, connect both graphs
    Add an edge between the nodes of the two graphs if they have:
    - colour graph in two colours as "missing" and "specific"
    - more than 2 descs in common
    - they are from the same cluster
    '''
    print("Connecting missing and specific graphs...")
    nx.set_node_attributes(specific, "specific", "type")
    nx.set_node_attributes(specific, "true", "property")
    nx.set_node_attributes(depleted, "missing", "type")
    nx.set_node_attributes(depleted, "true", "property")
    merged_graph = nx.union(specific, depleted, rename=("S-", "M-"))

    specific_type = {} # keep track of the cluster specific genes to see if they are fake, truncated, longer, same function or completely specific
    ## from the specific can infer for missing...

    for n1 in specific.nodes(data=True):
        specific_type[n1[0]] = {"missing": {}, "type":"truly cluster specific", "desc":"-"}
        for n2 in depleted.nodes(data=True):
            if n1[1]['cluster'] != n2[1]['cluster']: ## don't add an edge between genes from different clusters
                continue
            n1_descs = n1[1]['desc'].split("/")
            n2_descs = n2[1]['desc'].split("/")
            intersection = list(set(n1_descs) & set(n2_descs))
            if len(intersection) > 1: ## add an edge if they share some function
                merged_graph.add_edge("S-"+n1[0], "M-"+n2[0], weight = len(intersection), terms = "/".join(intersection))
                merged_graph = set_specific_type(merged_graph, specific_type, n1[0], "same functional group", missing = n2[0])

    ## add wrong edges using blast output
    print("Looking for wrong/truncated/fused genes from blast...")
    with open("specific_against_db.tab") as f:
        for line in f:
            toks = line.strip().split()
            toks[0] = re.sub(r"[^a-zA-Z0-9_\*-]+", '', toks[0])
            toks[1] = re.sub(r"[^a-zA-Z0-9_\*-]+", '', toks[1])
            curr_m = "M-" + toks[1]
            curr_s = "S-" + toks[0]
            flag = False
            ## if they are from different clusters they can't be wrong
            # if merged_graph.nodes[curr_m]["cluster"] != merged_graph.nodes[curr_s]["cluster"]:
            #     continue

            if float(toks[2]) > 90 and int(toks[4]) - 20 <= int(toks[5]) and int(toks[5]) <= int(toks[4]) + 20:
                if not merged_graph.has_edge(curr_s, curr_m):
                    merged_graph.add_edge(curr_s, curr_m, weight = "1", terms = "BLAST")
                merged_graph = set_specific_type(merged_graph,specific_type, toks[0], "wrong", missing = toks[1])
            elif float(toks[2]) > 95 and max(int(toks[3])/float(toks[4]),int(toks[3])/float(toks[5])) >= 0.7:
                if not merged_graph.has_edge(curr_s, curr_m):
                    merged_graph.add_edge(curr_s, curr_m, weight = "1", terms = "BLAST-truncated")
                if int(toks[4]) < int(toks[5]):
                    merged_graph = set_specific_type(merged_graph,specific_type, toks[0], "truncated", missing = toks[1])
                else:
                    merged_graph = set_specific_type(merged_graph,specific_type, toks[0], "longer", missing = toks[1])

    ccs = sorted(nx.connected_components(merged_graph), key=len, reverse=True)
    edges_to_remove = [] ## remove edges when there are only a small number of shared terms
    index = 0
    print("Splitting clusters using dbscan")
    for component in ccs:
        if len(component) < 2:
            continue
        ## remove terms which are rare in the connected component...
        terms = {}
        num_edges = 0
        for edge in merged_graph.edges(component, data = True):
            num_edges += 1
            curr_terms = edge[2]["terms"].split("/")
            for t in curr_terms:
                if t not in terms:
                    terms[t] = 0
                terms[t] += 1
        limit = np.quantile(np.array(terms.values()), 0.7)
        ## run dbscan on the cluster to seperate into groups that share a lot of terms
        ## remove edges
        curr_edges_to_remove, merged_graph = split_cluster(nx.Graph(merged_graph.subgraph(component)), merged_graph)
        edges_to_remove += curr_edges_to_remove


    nx.write_gml(merged_graph, "pre_merged_graph.gml")
    print("Pre merged graph written and complete...")
    merged_graph.remove_edges_from(edges_to_remove)
    ## create the details files:
    with open("specific_gene_types.csv", "w") as out1:
        with open("missing_gene_types.csv","w") as out2:
            out1.write("Gene,Type,Partner,Desc\n")
            out2.write("Gene,Type,Partner,Desc\n")
            for node in merged_graph.nodes(data = True):
                if node[0].startswith("M"):
                    out = out2
                    partner = "NA"
                else:
                    out = out1
                    if node[1]["property"] not in specific_type[node[0][2:]]["missing"]:
                        partner = "NA"
                    else:
                        partner = "/".join(list(specific_type[node[0][2:]]["missing"][node[1]["property"]]))

                out.write(node[0][2:] + "," + node[1]["property"] + "," + partner + "," + node[1]["desc"].replace(",","-") + "\n")
    ## create the functional clusters output, but remove
    ## nodes that are wrong or just truncated versions (pseudogenised? wrong chosen?)
    out = open("functional_clusters.csv","w")
    out.write("CC,Size,Num_clusters,Clusters,Type,Common_terms\n")
    remove = []
    for n in merged_graph.nodes(data=True):
        if n[1]["property"] in ["wrong", "truncated","longer"]:
            remove.append(n[0])
    merged_graph.remove_nodes_from(remove)
    nx.write_gml(merged_graph, "post_merged_graph.gml")
    return


def get_seqs():
    seqs = {}
    with open("db_for_specific_genes.fa") as handle:
        for values in SimpleFastaParser(handle):
            seqs[re.sub(r"[^a-zA-Z0-9_\*-]+", '', values[0].split()[0])] = values[1]
    with open("specific.fa") as handle:
        for values in SimpleFastaParser(handle):
            seqs[re.sub(r"[^a-zA-Z0-9_\*-]+", '', values[0].split()[0])] = values[1]
    return seqs

def generate_msa_output(filein):
    ''' create multiple sequence alignments for all genes
    that have some sort of relatioship'''
    seqs = get_seqs()
    with open(filein) as f:
        for line in f:
            if line.startswith("Gene"):
                continue
            toks = line.strip().split(",")
            if toks[2] == "NA":
                continue
            name = toks[0]
            partners = toks[2].split("/")
            fasta_out  = open(os.path.join("groups/", name + ".fa"), "w")
            fasta_out.write(">" + name + "(specific)\n" + seqs[name] + "\n")
            for p in partners:
                fasta_out.write(">" + name + "(missing)\n" + seqs[p] + "\n")
            fasta_out.close()
            out_msa = open(os.path.join("groups", name + "_msa.fa"), "w")
            p = subprocess.Popen(["mafft", os.path.join("groups", name + ".fa")], stdout=out_msa, stderr=subprocess.PIPE)
            p.wait()
            out_msa.close()
    return

if __name__ == "__main__":

    gene_classification, header = get_gene_classification()

    all_specific_terms = get_detailed_functions("specific_genes.csv", "specific_genes_detailed.csv")
    G1 = summarise_terms_as_network(all_specific_terms, "specific_genes_graph.gml")

    all_depleted_terms = get_detailed_functions("missing_genes.csv", "missing_genes_detailed.csv")
    G2 = summarise_terms_as_network(all_depleted_terms, "missing_genes_graph.gml")

    connect_both_graphs(G1, G2)
    generate_msa_output("specific_gene_types.csv")
    quit()
