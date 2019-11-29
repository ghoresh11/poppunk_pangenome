from Bio.SeqIO.FastaIO import SimpleFastaParser
from numpy import mean,std
import os

def read_gff_file(gff_file, member_props):
    ''' read a gff file and for each member return the contig length and
    distance from the edge'''
    contig_lengths = {}
    curr_contigs = {}
    num_rare_variants = 0

    with open(gff_file) as f:
        for line in f:
            if line.startswith("##sequence-region"):
                toks = line.strip().split()
                contig_lengths[toks[1]] = int(toks[3])
                continue
            if line.startswith("##FASTA"):
                break
            if line.startswith("#"):
                continue
            toks = line.strip().split("\t")
            if toks[2] != "CDS":
                continue
            gene_name = toks[-1].split(";")[0].replace("ID=", "")
            contig_length = contig_lengths[toks[0]]
            distance = min(contig_length - float(toks[4]), float(toks[3]))
            member_props[gene_name] = {"contig_length": contig_length,
                                       "distance": distance}
    return

member_props = {}
# step 1: get a list of all the GFF files
gff_files = []
with open("/lustre/scratch118/infgen/team216/gh11/e_coli_collections/FILTERED_MD_FINAL_ALL.tab") as f:
    for line in f:
        toks = line.strip().split("\t")
        if line.startswith("ID"):
            annot_loc = toks.index("New_annot_loc")
            continue
        read_gff_file(toks[annot_loc], member_props)

out = open("harder_gene_props.csv", "w")
out.write("Gene,Property,Mean,Std\n")
gene_cluster_loc = "/lustre/scratch118/infgen/team216/gh11/231019_clusters/"
all_genes = os.listdir(gene_cluster_loc)
for gene_file in all_genes:
    if not gene_file.endswith(".fa"):
        continue
    gene = gene_file.split(".")[0]
    gene_props = {}
    for key in ["contig_length","distance"]:
        gene_props[key] = []
    curr_fasta_file = os.path.join(gene_cluster_loc, gene_file)
    with open(curr_fasta_file) as handle:
        for values in SimpleFastaParser(handle):
            curr_member = values[0]
            if values[0] not in member_props: ## shouldn't happen in general
                continue
            for key in ["contig_length","distance"]:
                gene_props[key].append(member_props[curr_member][key])
    for prop in gene_props:
        out.write(",".join([gene, prop, str(mean(gene_props[prop])), str(std(gene_props[prop]))]) + "\n")

out.close()
