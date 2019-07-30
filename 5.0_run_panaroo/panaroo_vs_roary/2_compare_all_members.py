import os
import time


def clusterid_to_annotid(gene_data_file):
    ''' the panaroo ids are different from the annotation, get the diff'''
    print("Getting ids for %s..." %gene_data_file)
    cluster_to_annot = {}
    with open(gene_data_file) as f:
        for line in f:
            toks = line.strip().split(",")
            if line.startswith("gff_file"):
                cluster_id_loc = toks.index("clustering_id")
                annot_id_loc = toks.index("annotation_id")
                continue
            cluster_to_annot[toks[cluster_id_loc]] = toks[annot_id_loc]
    return cluster_to_annot


def init_network(roary_clustered_proteins_file):
    ''' create a graph with edges between members based on roary outputs'''
    print("Initiating graph according to roary outputs...")
    edges = {}
    with open(roary_clustered_proteins_file) as f:
        for line in f:
            toks = line.strip().split()[1:]
            for i in range(0, len(toks)-1):
                for j in range(i+1, len(toks)):
                    edge = [toks[i], toks[j]]
                    edge.sort()
                    edge = ";".join(edge)
                    edges[edge] = "roary"
    return edges

def add_panaroo_to_network(edges, ids, panaroo_presence_absence_file):
    ''' add the panaroo results to the graph'''
    print("Adding panaroo results to graph...")
    with open(panaroo_presence_absence_file) as f:
        for line in f:
            toks = line.strip().split(",")
            if line.startswith("Gene"):
                genome_indexes = toks.index("Avg group size nuc") + 1
                continue
            curr_members = toks[genome_indexes:]
            for i in range(0, len(curr_members) - 1):
                for j in range(i+1, len(curr_members)):
                    if curr_members[i] == "" or curr_members[j] == "":
                        continue
                    if "refound" in curr_members[i]:
                        member1 = curr_members[i]
                    else:
                        member1 = ids[curr_members[i]]
                    if "refound" in curr_members[j]:
                        member2 = curr_members[j]
                    else:
                        member2 = ids[curr_members[j]]

                    edge = [member1, member2]
                    edge.sort()
                    edge = ";".join(edge)

                    if edge in edges:
                        edges[edge] = "both"
                    else:
                        edges[edge] = "panaroo"
    return

roary_dir = "/lustre/scratch118/infgen/team216/gh11/e_coli_collections/poppunk/new_roary/"
panaroo_dir = "/lustre/scratch118/infgen/team216/gh11/e_coli_collections/poppunk/panaroo/"

out = open("summary_all.csv", "w")
out.write("cluster, both, roary_only, panaroo_only, jaccard\n")

for i in range(5,6):
    if i == 50:
        continue
    edges = init_network(os.path.join(roary_dir, str(i), "clustered_proteins"))
    ids = clusterid_to_annotid(os.path.join(panaroo_dir, str(i), "gene_data.csv"))
    add_panaroo_to_network(edges, ids, os.path.join(panaroo_dir, str(i), "gene_presence_absence.csv"))
    ## edges that are 00 or 11 are in both, 01 and 10 are in one but not in the other!
    ## 00 and 11 are the number of edges that are the same
    print("Summing values")
    both = sum(value == "both" for value in edges.values())
    roary = sum(value == "roary" for value in edges.values())
    panaroo = sum(value == "panaroo" for value in edges.values())
    jaccard = float(both) / (both + roary + panaroo)
    out.write(",".join(map(str, [i, both, roary, panaroo, jaccard])) + "\n")
    with open(str(i) + "_edge_list.txt", "w") as out2:
        for e in edges:
            out2.write("\t".join(e.split(";") + [edges[e]]) + "\n")

out.close()
