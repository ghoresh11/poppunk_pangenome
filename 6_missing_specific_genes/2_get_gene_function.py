import os
import re
import networkx as nx
from Bio.SeqIO.FastaIO import SimpleFastaParser
import subprocess


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
                            "Amino acid"]:
                                sigs.append(s)
                    all_gene_domains_sigs[g] = {"sigs":sigs, "cluster":toks[0]}
    return all_gene_domains_sigs


def summarise_terms_as_network(all_gene_terms, out):
    '''create a network with edges between the nodes if they share
    any functional terms'''
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


def get_seqs():
    seqs = {}
    with open("missing.fa") as handle:
        for values in SimpleFastaParser(handle):
            seqs[re.sub(r"[^a-zA-Z0-9_\*-]+", '', values[0].split()[0])] = values[1]
    with open("specific.fa") as handle:
        for values in SimpleFastaParser(handle):
            seqs[re.sub(r"[^a-zA-Z0-9_\*-]+", '', values[0].split()[0])] = values[1]
    return seqs

def set_specific_type(merged_graph, specific_type, key, new_type, missing = "NA",desc = "-"):
    ''' this is complicated, but there's a hierarchy in defining the property of each
    missing or cluster specific gene, depending on their relationships
    Basically wrong>truncated/missing>same functional group'''
    s_key = "S-" + key
    m_key = "M-" + missing
    orig = specific_type[key]["type"]
    if new_type == "wrong" or specific_type[key]["type"] == "wrong":
        specific_type[key]["type"] = "wrong"
        merged_graph.nodes[s_key]["property"] = "wrong"
        merged_graph.nodes[m_key]["property"] = "wrong"
        if orig != specific_type[key]["type"]:
            specific_type[key]["missing"] = missing
        return merged_graph

    if new_type in ["longer","truncated"]:
        specific_type[key]["type"] = new_type

    if specific_type[key]["type"] in ["longer","truncated"]:
        merged_graph.nodes[s_key]["property"] = specific_type[key]["type"]

        if merged_graph.nodes[m_key]["property"] != "wrong":
            if specific_type[key]["type"] == "longer":
                merged_graph.nodes[m_key]["property"] = "truncated"
            else:
                merged_graph.nodes[m_key]["property"] = "longer"
        if orig != specific_type[key]["type"]:
            specific_type[key]["missing"] = missing
        return merged_graph
    if new_type == "same functional group" or specific_type[key]["type"] == "same functional group":
        specific_type[key]["type"] = "same functional group"
        specific_type[key]["desc"] = desc
        merged_graph.nodes[s_key]["property"] =  "same functional group"
        if merged_graph.nodes[m_key]["property"] not in ["longer","truncated","wrong"]:
            merged_graph.nodes[m_key]["property"] =  "same functional group"
    if orig != specific_type[key]["type"]:
        specific_type[key]["missing"] = missing
    return merged_graph

def connect_both_graphs(specific, depleted, outdir):
    ''' in order to understand whether there are genes in a cluster that
    are cluster specific because they are compensated by a different gene with the
    same function, connect both graphs
    Add an edge between the nodes of the two graphs if they have:
    - colour graph in two colours as "missing" and "specific"
    - more than 2 descs in common
    - they are from the same cluster
    '''
    nx.set_node_attributes(specific, "specific", "type")
    nx.set_node_attributes(specific, "true", "property")
    nx.set_node_attributes(depleted, "missing", "type")
    nx.set_node_attributes(depleted, "true", "property")
    merged_graph = nx.union(specific, depleted, rename=("S-", "M-"))

    specific_type = {} # keep track of the cluster specific genes to see if they are fake, truncated, longer, same function or completely specific
    ## from the specific can infer for missing...

    for n1 in specific.nodes(data=True):
        specific_type[n1[0]] = {"missing":"NA","type":"truly cluster specific", "desc":"-"}
        for n2 in depleted.nodes(data=True):
            if n1[1]['cluster'] != n2[1]['cluster']:
                continue
            n1_descs = n1[1]['desc'].split("/")
            n2_descs = n2[1]['desc'].split("/")
            intersection = list(set(n1_descs) & set(n2_descs))
            if len(intersection) > 1:
                merged_graph.add_edge("S-"+n1[0], "M-"+n2[0], weight = len(intersection), terms = "/".join(intersection))

    if not os.path.exists(outdir):
        os.makedirs(outdir)
    ccs = sorted(nx.connected_components(merged_graph), key=len, reverse=True)
    seqs = get_seqs()
    edges_to_remove = [] ## remove edges when there are only a small number of shared terms
    index = 0
    for component in ccs:
        if len(component) < 2:
            continue
        ## get a name for this cluster
        names = []
        for node in component:
            names.append(node[2:6])
        name = max(set(names), key=names.count) + "_" + str(index)
        index += 1
        fasta_out  = open(os.path.join(outdir, name + ".fa"), "w")
        ## get members of each PopPUNK cluster
        for node in component:
            curr_cluster = merged_graph.nodes[node]["cluster"]
            fasta_out.write(">" + node + " (" + curr_cluster + ")" + "\n" + seqs[node[2:]] + "\n")
        fasta_out.close()
        out_msa = open(os.path.join(outdir, name + "_msa.fa"), "w")
        p = subprocess.Popen(["mafft", os.path.join(outdir, name + ".fa")], stdout=out_msa, stderr=subprocess.PIPE)
        p.wait()
        out_msa.close()

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
        terms_to_remove = []
        for t in terms:
            terms[t] = terms[t]/float(num_edges)
            if terms[t] < 0.8: ## this is a rare term in this cluster
                terms_to_remove.append(t)
        for edge in merged_graph.edges(component, data = True):
            curr_terms = edge[2]["terms"].split("/")
            final_terms = [x for x in curr_terms if x not in terms_to_remove]
            if len(final_terms) < 1: ## if remove random terms, the edge should be removed
                edges_to_remove.append((edge[0], edge[1]))
            else:
                edge[2]["terms"] = "/".join(final_terms) ## remove some of the terms
        ## go over the MSA file and see which genes are:
        # 1. mistakes, almost identical and were split for some reason
        # 2. Truncated, which could also be mistakes because of how roary chooses the rep
        # 3. Truly different
        aligned_seqs = {}
        with open(os.path.join(outdir, name + "_msa.fa")) as handle:
            for values in SimpleFastaParser(handle):
                aligned_seqs[values[0]] = values[1]
        ordered_seqs = list(aligned_seqs.keys())

        for i in range(0, len(ordered_seqs)-1):
            for j in range(i, len(ordered_seqs)):
                if ordered_seqs[i][0] == ordered_seqs[j][0]: # both specific or both missing
                    continue

                ## skip if they are from different clusters, they don't count in the relationship!
                if ordered_seqs[i].split()[-1] != ordered_seqs[j].split()[-1]:
                    continue


                seqa = aligned_seqs[ordered_seqs[i]]
                seqb = aligned_seqs[ordered_seqs[j]]
                gene_i = ordered_seqs[i][2:].split()[0]
                gene_j = ordered_seqs[j][2:].split()[0]
                node_i = ordered_seqs[i].split()[0]
                node_j = ordered_seqs[j].split()[0]

                if not merged_graph.has_edge(node_i, node_j): ## doesn't make sense they should have a relationship if dont share function (=edge)
                    continue

                diff = float(sum(1 for a, b in zip(seqa, seqb) if a != b)) ## num SNPs difference
                gaps = float(sum(1 for a, b in zip(seqa, seqb) if (a != b) and (a == "-" or b == "-")))

                len_seqa = len(seqa.replace("-",""))
                len_seqb = len(seqb.replace("-",""))

                if gaps > 0 and gaps >= (diff*0.95): ## almost all differences are due to gaps
                    if ordered_seqs[i][0] == "S":
                        if len_seqa > len_seqb:
                            merged_graph = set_specific_type(merged_graph,specific_type, gene_i, "longer", missing = gene_j)
                        else:
                            merged_graph = set_specific_type(merged_graph,specific_type, gene_i, "truncated", missing = gene_j)
                    else:
                        if len_seqb > len_seqa:
                            merged_graph = set_specific_type(merged_graph,specific_type, gene_j, "longer", missing = gene_i)
                        else:
                            merged_graph = set_specific_type(merged_graph,specific_type, gene_j, "truncated", missing = gene_i)
                elif diff/max(len_seqa, len_seqb) <= 0.1: ## they're 90% similar, this is an error
                    if ordered_seqs[i][0] == "S":
                        merged_graph = set_specific_type(merged_graph,specific_type, gene_i, "wrong", missing = gene_j)
                    else:
                        merged_graph = set_specific_type(merged_graph,specific_type, gene_j, "wrong", missing = gene_i)
                else: ## the two genes have functional relationship
                    ## only change it if it hasn't been changed, there's a hierarchy of labels wrong > truncated/longer > same functional category
                    if ordered_seqs[i][0] == "S":
                        merged_graph = set_specific_type(merged_graph,specific_type, gene_i, "same functional group", missing = gene_j, desc = merged_graph[ordered_seqs[i].split()[0]][ordered_seqs[j].split()[0]]["terms"].replace(",","/"))
                    else:
                        merged_graph = set_specific_type(merged_graph,specific_type, gene_j, "same functional group", missing = gene_i, desc = merged_graph[ordered_seqs[i].split()[0]][ordered_seqs[j].split()[0]]["terms"].replace(",","/"))
    nx.write_gml(merged_graph, "pre_merged_graph.gml")
    merged_graph.remove_edges_from(edges_to_remove)
    with open("specific_gene_types.csv","w") as out:
        out.write("Specific_gene,Missing_gene,Type,Desc\n")
        for gene in specific_type:
            out.write(",".join([gene, specific_type[gene]["missing"], specific_type[gene]["type"], specific_type[gene]["desc"]]) + "\n")
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
    index = 0
    if not os.path.exists(outdir + "_2"):
        os.makedirs(outdir+ "_2")
    for component in sorted(nx.connected_components(merged_graph), key=len, reverse=True):
        fasta_out  = open(os.path.join(outdir, name + ".fa"), "w")
        ## get members of each PopPUNK cluster
        for node in component:
            curr_cluster = merged_graph.nodes[node]["cluster"]
            fasta_out.write(">" + node + " (" + curr_cluster + ")" + "\n" + seqs[node[2:]] + "\n")
        fasta_out.close()
        out_msa = open(os.path.join(outdir, name + "_msa.fa"), "w")
        p = subprocess.Popen(["mafft", os.path.join(outdir, name + ".fa")], stdout=out_msa, stderr=subprocess.PIPE)
        p.wait()
        out_msa.close()

        size = len(component)
        clusters = set()
        types = set()
        for node in component:
            clusters.add(merged_graph.nodes[node]["cluster"])
            types.add(merged_graph.nodes[node]["type"])
        clusters = list(clusters)
        num_clusters = len(clusters)
        if len(types) == 1:
            types = list(types)[0]
        else:
            types = "both"
        ## get the terms in this cluster
        terms = {}
        for edge in merged_graph.edges(component, data = True):
            curr_terms = edge[2]["terms"].split("/")
            for t in curr_terms:
                if t not in terms:
                    terms[t] = 0
                terms[t] += 1
        terms_out = []
        for t in terms:
            terms_out.append(t + ":" + str(terms[t]))
        index += 1
        out.write(",".join(map(str, [index, size, num_clusters, "/".join(clusters),types,  "/".join(terms_out).replace(",","-").replace("\n","")])) + "\n")
    out.close()
    return merged_graph

def create_heatmap_file(merged_graph):
    ''' for each of the functional groups, generate a heatmap with
    presence and absence of that functional group with in a lineage
    the file should have the functions as rows, the clusters as columns
    and should have a count for each cluster in each functional group'''
    out = open("heatmap.csv","w")
    clusters_out = []
    for i in range(1,52):
        if i not in [21,43,49]:
            clusters_out.append(i)
    out.write("Group," + ",".join(map(str, clusters_out)) + "\n")
    for component in sorted(nx.connected_components(merged_graph), key=len, reverse=True):
        curr = [0] * len(clusters_out)
        terms = set()
        for edge in merged_graph.edges(component, data = True):
            curr_terms = edge[2]["terms"].split("/")
            for t in curr_terms:
                terms.add(t)
        for node in component:
            curr_cluster = merged_graph.nodes[node]["cluster"]
            if  merged_graph.nodes[node]["type"] == "specific":
                curr[clusters_out.index(curr_cluster)] += 1
        out.write("/".join(terms) + "," + ",".join(map(str,curr)) + "\n")
    out.close()
    return

if __name__ == "__main__":
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

    # #get_detailed_functions("depleted_genes.csv", "depleted_genes_detailed.csv")
    all_specific_terms = get_detailed_functions("specific_genes.csv", "specific_genes_detailed.csv")
    G1 = summarise_terms_as_network(all_specific_terms, "specific_genes_graph.gml")

    all_depleted_terms = get_detailed_functions("missing_genes.csv", "missing_genes_detailed.csv")
    G2 = summarise_terms_as_network(all_depleted_terms, "missing_genes_graph.gml")
    merged_graph = connect_both_graphs(G1, G2, "merged")
    create_heatmap_file(merged_graph)
