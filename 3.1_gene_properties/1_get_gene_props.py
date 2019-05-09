import argparse
import os
from numpy import mean,std
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.SeqUtils import GC

''' for each gene cluster in all clusters, get all the members of that cluster
and calculate the following properties:
1. Mean length + SD
2. Mean GC content + SD
3. Mean relative location on the chromosome + SD
4. Mean contig length + SD
These prorperties can later be used for the following:
1. Understanding the properties of the gene clusters that are later merged to be
defined as the same, i.e. are they usually short?
2. Look at the properties of core, inter and rare genes relative to each other
3. Look at the properties of genes from the different functional families.
'''

def get_roary_dir(args):
    ''' return the roary directory for this cluster'''
    directories = os.listdir(args.roary_dir)
    for d in directories:
        curr_cluster = d.split("_")[0]
        if curr_cluster == str(args.cluster):
            input_dir = os.path.join(args.roary_dir, d)
            break
    return input_dir

def get_genomes_to_gene_clusters(args):
    ''' read the clustered proteins output.
    Return a dictionary of genome -> name_of_gene_in_genome -> name of geneCluster
    '''
    print("Reading clustered proteins file...")
    input_dir = get_roary_dir(args)
    genomes_gene_cluster = {}
    with open(os.path.join(input_dir, "clustered_proteins")) as f:
        for line in f:
            toks = line.strip().split()
            gene_cluster = toks[0].replace(":","")
            for gene in toks[1:]:
                genome = "_".join(gene.split("_")[:-1])
                if genome not in genomes_gene_cluster:
                    genomes_gene_cluster[genome] = {}
                genomes_gene_cluster[genome][gene] = gene_cluster

    return genomes_gene_cluster


def get_props_for_one_genome(args, properties, genome_name, gene_cluster, gene_reps):
    ''' open the GFF file for a single genome, and obtain all the
    properties for the gene clusters of this file
    Return: nothing
    Do: udated the gene cluster properites based on this file'''
    ## get location of GFF file
    with open(os.path.join(args.gff_list_dir, "jobs_" + str(args.cluster) + ".txt")) as f:
        for line in f:
            if genome_name in line:
                gff_file = line.strip()
                break

    contig_lengths = {}
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
            if gene_name not in gene_cluster:
                ## some genes are not in the final roary output
                continue

            curr_cluster = gene_cluster[gene_name]
            if curr_cluster not in properties:
                properties[curr_cluster] = {"class": "NA",
                                            "cog" : "NA",
                                            "GC": "NA",
                                            "contig_lengths":[],
                                            "positions": [],
                                            "protein_length": []}

            if args.cluster == 1:
                gene_name = curr_cluster    # special case because the output is different

            if gene_name in gene_reps:
                properties[curr_cluster]["class"] = gene_reps[gene_name]["class"]
                properties[curr_cluster]["cog"] = gene_reps[gene_name]["cog"]
                properties[curr_cluster]["GC"] = gene_reps[gene_name]["GC"]

            contig_length = contig_lengths[toks[0]]
            properties[curr_cluster]["contig_lengths"].append(contig_length)
            if toks[6] == "+":
                end = float(toks[4])
            else:
                end = contig_length - float(toks[3]) ## minus strand
            properties[curr_cluster]["positions"].append(contig_length - end)
            properties[curr_cluster]["protein_length"].append((int(toks[4]) - int(toks[3]) + 1) / 3.0)
    return

def read_eggnog(args):
    ''' get the class and COG from the eggnog outputs for each gene cluster
    return the reps class and COG category'''
    print("Getting COG details...")
    all_files = os.listdir(args.eggnog_dir)
    files = {}
    for f in all_files:
        if f.split("_")[0] == str(args.cluster):
            gene_class = f.split(".")[0].replace(str(args.cluster) + "_", "")
            files[gene_class] = os.path.join(args.eggnog_dir, f)

    gene_reps = {}
    for gene_class in files:
        with open(files[gene_class]) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                toks = line.strip().split("\t")
                gene_name = toks[0]
                cog = "?"
                if len(toks) > 12:
                    cog = toks[-2].replace(",", "/").replace(" ", "")
                gene_reps[gene_name] = {"cog": cog, "class": gene_class}
    return gene_reps


def get_rep_class(args, gene_reps):
    ''' get the class of all the genes using the outputs in the roary directories'''
    print("Getting class details....")
    input_dir = get_roary_dir(args)
    for gene_class in ["rare", "inter", "soft_core", "core"]:
        with open(os.path.join(input_dir, gene_class + "_genes.fa")) as handle:
            for values in SimpleFastaParser(handle):
                name = values[0].split()[0]
                if name not in gene_reps.keys():
                    gene_reps[name] = {}
                    gene_reps[name]["class"] = "?"
                    gene_reps[name]["cog"] = "?"


def get_rep_GC_content(args, gene_reps):
    ''' get the GC content of the reprasentative sequence for each gene cluster
    update this in the gene_reps dictionary'''
    print("Calculating GC content...")
    input_dir = get_roary_dir(args)
    with open(os.path.join(input_dir, "pan_genome_reference.fa")) as handle:
        for values in SimpleFastaParser(handle):
            name = values[0].split()[0]
            if name in gene_reps.keys():
                gene_reps[name]["GC"] = GC(values[1])
            else:
                gene_reps[name] = {}
                gene_reps[name]["GC"] = GC(values[1])
                gene_reps[name]["class"] = "?"
                gene_reps[name]["cog"] = "?"
    return



def get_properties_for_all(args, genome_gene_cluster):
    ''' go over other outputs and GFF files to obtain the properties of all the genes '''

    ### For COG category and GC content, I'm not using an average but using the rep
    ## from the roary outputs. This is problematic with cluster 1 because I changed the
    # output for that cluster.
    gene_reps = read_eggnog(args)
    get_rep_class(args, gene_reps)
    get_rep_GC_content(args, gene_reps)

    properties = {}
    print("Calculating properties per genome...")
    for genome in genome_gene_cluster:
        get_props_for_one_genome(args, properties, genome, genome_gene_cluster[genome], gene_reps)

    print("Generating output...")
    with open(os.path.join(str(args.cluster) + "_gene_properties.csv"), "w") as out:
        out.write("gene, class, COG, cnt, GC,  mean_length, sd_length, mean_pos, sd_pos, mean_contig_length, sd_contig_length\n")
        for gene in properties:
            props = properties[gene]
            out.write(",".join(map(str,
            [gene, props["class"], props["cog"], len(props["protein_length"]), props["GC"],
            mean(props["protein_length"]),  std(props["protein_length"]),
            mean(props["positions"]),  std(props["positions"]), mean(props["contig_lengths"]),  std(props["contig_lengths"])])) + "\n")

    return


def run(args):
    genomes_gene_cluster = get_genomes_to_gene_clusters(args)
    get_properties_for_all(args, genomes_gene_cluster)
    print("DONE!")
    return


def get_options():
    parser = argparse.ArgumentParser(description='Extract the gene sequences from roary outputs, and merge rare genes')
    parser.add_argument('--cluster', required=False,
                        type=int, default = 30,
                        help='Cluster for which genome properties are being obtained [%(default)s]')
    parser.add_argument('--gff_list_dir', required=False,
                        type=str, default = "/lustre/scratch118/infgen/team216/gh11/e_coli_collections/poppunk/dists_analysis/gff_jobs",
                        help='Path to directory with all the GFF files [%(default)s]')
    parser.add_argument('--roary_dir', required=False,
                        type=str, default = "/lustre/scratch118/infgen/team216/gh11/e_coli_collections/poppunk/dists_analysis/roary_outputs",
                        help='Path to directory with all the original roary outputs [%(default)s]')
    parser.add_argument('--eggnog_dir', required=False,
                        type=str, default = "/lustre/scratch118/infgen/team216/gh11/e_coli_collections/poppunk/eggnog",
                        help='Path to directory with all the EGGNOG outputs [%(default)s]')
    return parser.parse_args()


if __name__ == "__main__":
    # get arguments from user
    options = get_options()
    run(options)


### to run
