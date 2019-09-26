import os
from Bio.SeqIO.FastaIO import SimpleFastaParser
import re
from Bio.Seq import translate

''' extract all the gene sequences from the original roary outputs
and save a fasta file for each gene class to be analysed for function
and to find which genes are missing or enriched in specific clusters'''


## step 1-> get the cluster members of each gene, better to point from cluster -> member -> to gene
print('Reading members file....')
cluster_member_gene = {}
with open("/lustre/scratch118/infgen/team216/gh11/e_coli_collections/poppunk/new_roary/pairwise/analysis/190919/members.csv") as f:
    for line in f:
        if line.startswith("Gene"):
            continue
        toks = line.strip().split(",")
        gene = toks[0]
        members = toks[1].split("\t")
        for m in members:
            curr_cluster = m.split("_")[0]
            curr_member = "_".join(m.split("_")[1:])
            if curr_cluster not in cluster_member_gene:
                cluster_member_gene[curr_cluster] = {}
            cluster_member_gene[curr_cluster][curr_member] = gene

## step 2 -> go over all the pan_reference_genome files and choose the longest sequence for each gene
print("Getting sequenced from pan-genome-reference files...")
gene_sequence = {}
for cluster in cluster_member_gene:
    print(cluster)
    curr_members_gene = cluster_member_gene[cluster]
    curr_ref_file = os.path.join("/lustre/scratch118/infgen/team216/gh11/e_coli_collections/poppunk/new_roary", cluster, "longest_pan_genome_reference.fa")
    with open(curr_ref_file) as handle:
        for values in SimpleFastaParser(handle):
            member = values[0].split()[-1]
            curr_gene = curr_members_gene[member]
            if curr_gene not in gene_sequence or len(gene_sequence[curr_gene]["seq"]) < len(values[1]):
                gene_sequence[curr_gene] = {"seq": values[1], "rep": cluster + "_" + member}

## step 3 -> read the gene classification file and create a file for each class, can concat all later.
print("Generating outputs...")
files_out = {}
#with open("gene_classification.csv") as f:
with open("soft_core_genes.txt") as f:
    for line in f:
        toks = line.strip().split("\t")
        desc = "190919_soft_core_genes"
        if len(toks)>1:
            desc = toks[1]
        as_file_name = re.sub('[^a-zA-Z0-9\n\._]', '_', desc)
        if as_file_name not in files_out:
            files_out[as_file_name] = open(as_file_name + ".fa", "w")
            files_out[as_file_name + "_prot"] = open(as_file_name + "_prot.fa", "w")
        files_out[as_file_name].write(">" + toks[0] + " " + gene_sequence[toks[0]]["rep"] +  "\n" + gene_sequence[toks[0]]["seq"] + "\n")
        files_out[as_file_name + "_prot"].write(">" + toks[0] + " " + gene_sequence[toks[0]]["rep"] +  "\n" + translate(gene_sequence[toks[0]]["seq"]) + "\n")
for file_out in files_out:
    files_out[file_out].close()


## run interpro-scan
# for f in *.fa;
# do farm_interproscan -a $f -o ${f}_result;
# done
