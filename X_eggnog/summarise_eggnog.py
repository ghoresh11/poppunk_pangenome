import os
from Bio.SeqIO.FastaIO import SimpleFastaParser

def get_input_dirs(input_dir):
    ''' check all directories in the input dir
    return a list of directories that have all the required files'''
    print("Getting the input directories...")
    input_dir = os.path.abspath(input_dir)
    directories = [x[0] for x in os.walk(input_dir)]
    dirs_to_return = []
    for d in directories:
        if os.path.isfile(os.path.join(d, "pan_genome_reference.fa")) and os.path.isfile(os.path.join(d, "gene_presence_absence.csv")) and os.path.isfile(os.path.join(d, "gene_presence_absence.Rtab")):
            # for debugging, uncomment:
            # if d.split("/")[-1].split("_")[0] != "39":
            #      continue
            dirs_to_return.append(d)
    return dirs_to_return


def get_genes(input_dirs, gene_types):
    '''read all the genes from the .fa files into
    dictionary
    cluster -> -> type -> gene -> {"length": gene_length}'''
    print("Getting genes from the FASTA files...")
    genes = {}
    for d in input_dirs:
        cluster = d.split("/")[-1].split("_")[0]
        print(cluster)
        genes[cluster] = {}
        for gene_type in gene_types:
            genes[cluster][gene_type] = {}
            with open(os.path.join(d, gene_type + "_genes.fa")) as handle:
                for values in SimpleFastaParser(handle):
                    genes[cluster][gene_type][values[0].split()[0]] = {"length": len(values[1])}
    return genes

def read_eggnog_annotation(genes, gene_types, eggnog_dir):
    '''read the eggnog annotations and update
    the genes dictionary:
    cluster -> gene -> {"COG_ID":..., and whatever else}'''
    print("Reading the eggnot outputs..")
    for cluster in genes:
        print(cluster)
        for gene_type in gene_types:
            with open(os.path.join(eggnog_dir, cluster + "_" + gene_type + ".emapper.annotations")) as f:
                for line in f:
                    if line.startswith("#"):
                        continue
                    toks = line.strip().split("\t")
                    curr_dict = genes[cluster][gene_type][toks[0]]
                    if len(toks) >= 12:
                        curr_dict["COG"] = toks[-2].replace(",","/").replace(" ","")
    return


def generate_output(genes, gene_types):
    '''generate output for all (that can be easily plotted in R)
    cluster, gene_type, category, count'''
    print("Generating the output file...")
    out = open("eggnog_cog_summary.csv", "w")
    out.write("cluster, gene_type, cog_cat, count\n")
    for c in genes:
        for gene_type in gene_types:
            curr_counts = {}
            for gene in genes[c][gene_type]:
                if "COG" not in genes[c][gene_type][gene]:
                     genes[c][gene_type][gene]["COG"] = "?"
                curr_cogs = genes[c][gene_type][gene]["COG"].split("/")
                for cog in curr_cogs:
                    if cog not in curr_counts:
                        curr_counts[cog] = 0
                    curr_counts[cog] += 1
            for cog in curr_counts:
                out.write(",".join([c, gene_type, cog, str(curr_counts[cog])]) + "\n")
    out.close()
    return

def run():
    gene_types = ["rare", "inter", "core", "soft_core"]
    input_dirs = get_input_dirs("dists_analysis/roary_outputs/")
    genes = get_genes(input_dirs, gene_types)
    read_eggnog_annotation(genes, gene_types, "eggnog/")
    generate_output(genes, gene_types)
    return

if __name__ == "__main__":
    run()
