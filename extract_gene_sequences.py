import argparse
import os
from Bio.SeqIO.FastaIO import SimpleFastaParser
import subprocess
from Bio.Seq import translate


def get_input_dirs(input_dir):
    ''' check all directories in the input dir
    return a list of directories that have all the required files'''
    print("Getting the input directories...")
    input_dir = os.path.abspath(input_dir)
    directories = [x[0] for x in os.walk(input_dir)]
    dirs_to_return = []
    for d in directories:
        if os.path.isfile(os.path.join(d, "pan_genome_reference.fa")) and os.path.isfile(os.path.join(d, "gene_presence_absence.csv")) and os.path.isfile(os.path.join(d, "gene_presence_absence.Rtab")):
            dirs_to_return.append(d)
    return dirs_to_return


def classify_genes(input_dirs):
    ''' go over the Rtab file for each roary cluster,
    calculate the frequency of each gene in the cluster
    return: dict classification of genes into core, soft_core, intermediate, rare'''
    print("Calculating gene frequencies....")
    gene_class = {}
    for d in input_dirs:
        print(d)
        gene_class[d] = {}
        outputs = {}
        with open(os.path.join(d, "gene_presence_absence.Rtab")) as f:
            for line in f:
                if line.startswith("Gene"):
                    continue
                toks = line.strip().split("\t")
                freq = sum(map(int, toks[1:])) / float(len(toks)-1)
                gene_class[d][toks[0]] = freq
    return gene_class

def get_gene_sequences(input_dirs, gene_freqs):
    ''' use the dictionary, and the pan_genome_reference fasta file
    to create four different files for each cluster with the
    sequences of the genes in each category
    These will be used to postprocess the rare genes, and later
    to compare between different clusters.'''
    print("Getting the gene sequences....")
    out = open("Gene_lengths.csv", "w")
    out.write("Cluster, Type, Length\n")
    for d in input_dirs:
        cluster = d.split("/")[-1].split("_")[0]
        print(d)
        curr_cluster = gene_freqs[d]
        outputs = {}
        for type in ["rare", "inter", "soft_core", "core"]:
            outputs[type] = open(os.path.join(d, type + "_genes.fa"), "w")
        with open(os.path.join(d, "pan_genome_reference.fa")) as handle:
            for values in SimpleFastaParser(handle):
                gene_name = values[0].split()[-1]
                protein_sequence = translate(values[1])
                if curr_cluster[gene_name] < 0.15:
                    outputs["rare"].write(">" + values[0] + "\n" + protein_sequence + "\n")
                    out.write(cluster + ",rare, " + str(len(protein_sequence)) + "\n")
                elif curr_cluster[gene_name] < 0.95:
                    outputs["inter"].write(">" + values[0] + "\n" + protein_sequence + "\n")
                    out.write(cluster + ",inter, " + str(len(protein_sequence)) + "\n")
                elif curr_cluster[gene_name] < 0.99:
                    outputs["soft_core"].write(">" + values[0] + "\n" + protein_sequence + "\n")
                    out.write(cluster + ",soft_core, " + str(len(protein_sequence)) + "\n")
                else:
                    outputs["core"].write(">" + values[0] + "\n" + protein_sequence + "\n")
                    out.write(cluster + ",core, " + str(len(protein_sequence)) + "\n")
        for o in outputs:
            outputs[o].close()
    out.close()
    return

def summarise_functional_annotations(input_dirs, gene_freqs):
    ''' Use the presence absence CSV file to count the functional annotations
    of the genes from the different groups in each cluster'''
    print("Getting functional annotations....")
    for d in input_dirs:
        print(d)
        curr_cluster = gene_freqs[d]
        outputs = {}
        for type in ["rare", "inter", "soft_core", "core"]:
            outputs[type] = open(os.path.join(d, type + "_functions.txt"), "w")
        with open(os.path.join(d, "gene_presence_absence.csv")) as f:
            for line in f:
                toks = line.strip().split(",")
                if toks[0] == "\"Gene\"":
                    continue
                gene_name = toks[0].replace("\"","")
                annotation = toks[2].replace("\"","")
                if curr_cluster[gene_name] < 0.15:
                    outputs["rare"].write(annotation + "\n")
                elif curr_cluster[gene_name] < 0.95:
                    outputs["inter"].write(annotation + "\n")
                elif curr_cluster[gene_name] < 0.99:
                    outputs["soft_core"].write(annotation + "\n")
                else:
                    outputs["core"].write(annotation + "\n")
        for o in outputs:
            outputs[o].close()
    return


def filter_blast_results(input_file, output_file, min_identity):
    ''' filter the blast outputs to only include lines with matching min_identity'''
    with open(output_file, "w") as out:
        with open(input_file) as f:
            for line in f:
                ident = float(line.strip().split("\t")[-1])
                if ident > min_identity:
                    out.write(line)
    return


def blast_rare_genes(input_dirs, makeblastdb, blastp, cpus, mcl, inflation, min_identities):
    ''' Blast the rare genes of a file against each other
    return a dictionary with rare gene -> new rep'''
    print("Running blastp and MCL...")

    for d in input_dirs:
        print(d)

        subprocess.call([makeblastdb, "-in", os.path.join(d, "rare_genes.fa"), "-dbtype", "prot"])
        subprocess.call([blastp, "-query", os.path.join(d, "rare_genes.fa"),
                    "-db", os.path.join(d, "rare_genes.fa"), "-out", os.path.join(d, "rare_genes_blast_results.tab"), "-num_threads",
                    str(cpus), "-outfmt", "6 qseqid sseqid pident", "-evalue", "0.01"])


        for min_identity in min_identities:
            print("min identity: %d" %min_identity)
            new_reps = {}
            filtered_file = os.path.join(d, "rare_genes_blast_results_filtered_" + str(min_identity) + ".tab")
            mcl_file = os.path.join(d, "rare_genes_mcl_clusters_" + str(min_identity) + ".txt")
            filter_blast_results(os.path.join(d, "rare_genes_blast_results.tab"), filtered_file, min_identity)
            subprocess.call([mcl, filtered_file, "--abc", "-I", str(inflation) , "-o", mcl_file])
            with open(mcl_file) as f:
                for line in f:
                    toks = line.strip().split("\t")
                    rep = toks[0].split("|")[-1]
                    for i in range(0,len(toks)): ## keep the gene name, not the full rep
                        gene = toks[i].split("|")[-1]
                        new_reps[gene] = rep
            rewrite_pangenome_outputs(d, min_identity, new_reps)
    return new_reps


def rewrite_pangenome_outputs(d, min_identity, new_reps):
    ''' rewrite the roary output files so that the rare genes are merged into
    the clusters returned by new reps.
    The files that need to be changed are the "rare.fa" to include less genes
    The gene_presence_absence file also needs to merge rows which are all
    with the same rep now.'''
    ## 1. rewrite the FASTA output
    print("Rewriting output files...")
    out = open(os.path.join(d, "rare_genes_merged_" + str(min_identity) + ".fa"), "w")
    with open(os.path.join(d, "rare_genes.fa")) as handle:
        for values in SimpleFastaParser(handle):
            gene_name = values[0].split()[-1]
            if gene_name not in new_reps:
                print("PROBLEM WITH: %s (%s)!!!" %(gene_name, d)) ## will keep it
                out.write(">" + values[0] + "\n" + values[1] + "\n")
                continue
            if new_reps[gene_name] != gene_name: ## this gene is now represented by something else
                continue
            out.write(">" + values[0] + "\n" + values[1] + "\n")
    out.close()

    presence_absence = {}
    out = open(os.path.join(d, "gene_presence_absence_merged_" + str(min_identity) + ".Rtab"), "w")
    with open(os.path.join(d, "gene_presence_absence.Rtab")) as f:
        for line in f:
            if line.startswith("Gene"):
                out.write(line)
                continue
            toks = line.strip().split("\t")
            gene_name = toks[0]
            values = toks[1:]
            if gene_name not in new_reps: ## not a rare gene
                out.write(line)
                continue
            rep = new_reps[gene_name]
            if rep not in presence_absence:
                presence_absence[rep] = map(int, values)
            else:
                presence_absence[rep] = presence_absence[rep] + map(int, values) # add them
    for rep in presence_absence:
        presence_absence[rep] = [1 if x > 0 else 0 for x in presence_absence[rep]] # make 0 or 1
        out.write(rep + "\t" + "\t".join(map(str, presence_absence[rep])) + "\n")
    out.close()
    return

def summarise_outputs(input_dirs, min_identities):
    out = open("rare_filtering_summary.csv", "w")
    out.write("Cluster, Threshold, Gene_count, Decrease\n")

    for d in input_dirs:
        orig_gene_count = sum(1 for line in open(os.path.join(d,"rare_genes.fa"))) / 2 ## original number of genes
        cluster_num = d.split("/")[-1].split("_")[0]
        out.write(cluster_num + ",100," + str(orig_gene_count) + ",0\n")
        for min_identity in min_identities:
            merged_file = os.path.join(d, "rare_genes_merged_" + str(int(min_identity)) + ".fa")
            gene_count = sum(1 for line in open(merged_file)) / 2
            out.write(",".join(map(str,[cluster_num, min_identity, gene_count, (orig_gene_count - float(gene_count))/100.0 ])) + "\n")
    out.close()


def remove_temp_files(input_dirs, min_identities):
    for d in input_dirs:
        os.remove(os.path.join(d, "rare_genes.fa.phr"))
        os.remove(os.path.join(d, "rare_genes.fa.pin"))
        os.remove(os.path.join(d, "rare_genes.fa.psq"))
        # os.remove(os.path.join(d, "rare_genes_blast_results.tab"))
        for m in min_identities:
            ## remove all the temporary files, don't need them anymore
            os.remove(os.path.join(d, "rare_genes_blast_results_filtered_" + str(m) + ".tab"))
            os.remove(os.path.join(d, "rare_genes_mcl_clusters_"  + str(m) + ".txt"))
    return

def run(args):
    args.min_identities = map(float, args.min_identities.split(","))
    input_dirs = get_input_dirs(args.input_dir)
    gene_freqs = classify_genes(input_dirs)
    get_gene_sequences(input_dirs, gene_freqs)
    quit()
    summarise_functional_annotations(input_dirs, gene_freqs)
    blast_rare_genes(input_dirs, args.makeblastdb, args.blastp, args.cpu, args.mcl, args.inflation, args.min_identities)
    summarise_outputs(input_dirs, args.min_identities)
    remove_temp_files(input_dirs, min_identities)
    ## here: run eggnog to get the functions of each and create a summary for each cluster
    return


def get_options():
    parser = argparse.ArgumentParser(description='Extract the gene sequences from roary outputs, and merge rare genes')
    # input options
    parser.add_argument('--input_dir', required=False,
                        type=str, default =  ".",
                        help='path to input directory [%(default)s]')
    parser.add_argument('--blastp',
                        type=str, default = "/software/pubseq/bin/ncbi_blast+/blastp",
                        help='blastp executable [%(default)s]')
    parser.add_argument('--makeblastdb',
                        type=str, default = "/software/pubseq/bin/ncbi_blast+/makeblastdb",
                        help='makeblastdb executable [%(default)s]')
    parser.add_argument('--mcl',
                            type=str, default = "mcl",
                            help='mcl executable [%(default)s]')
    parser.add_argument('--cpu',
                        type=int, default = 16,
                        help='Number of CPUs to use [%(default)s]')
    parser.add_argument('--min_identities',
                            type=str, default = "50,60,70,80,90",
                            help='Comma seperated min identities to run with BLAST [%(default)s]')
    parser.add_argument('--inflation',
                                type=float, default = 1.5,
                                help='Inflation parameter to use in MCL [%(default)s]')
    return parser.parse_args()


if __name__ == "__main__":
    # get arguments from user
    options = get_options()
    run(options)
