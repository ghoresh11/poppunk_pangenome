import argparse
import os
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.Seq import reverse_complement
from Bio.Seq import translate


def merge_files(d):
    out_presence_absence = open(os.path.join(d, "merged_gene_presence_absence.Rtab"), "w")
    out_ref = open(os.path.join(d, "merged_pan_genome_reference.fa"), "w")
    out_cp = open(os.path.join(d, "merged_clustered_proteins"), "w")
    all_files = os.listdir(d)
    write_header = True
    for f in all_files:
        toks = f.split("_")
        ## only keep files of the format new0,1000
        if "new" not in toks[0] or toks[0] == "new" or toks[0] == "new2": # ignore old files
            continue
        skip_first = False
        if "clustered_proteins" in f:
            out = out_cp
        elif "pan_genome_reference" in f:
            out = out_ref
        else:
            out = out_presence_absence
            if not write_header:
                skip_first = True
            write_header = False
        with open(os.path.join(d,f)) as f_open:
            for line in f_open:
                if skip_first:
                    skip_first = False
                    continue
                out.write(line)
    out_presence_absence.close()
    out_ref.close()
    out_cp.close()
    return

def classify_genes(d):
    ''' go over the Rtab file for each roary cluster,
    calculate the frequency of each gene in the cluster
    return: dict classification of genes into core, soft_core, intermediate, rare'''
    print("Calculating gene frequencies....")
    gene_class = {}
    with open(os.path.join(d, "merged_gene_presence_absence.Rtab")) as f:
        for line in f:
            if line.startswith("Gene"):
                continue
            toks = line.strip().split("\t")
            freq = sum(map(int, toks[1:])) / float(len(toks)-1)
            gene_class[toks[0]] = freq
    output_genes_per_class(d, gene_class)
    return


def output_genes_per_class(d, gene_freqs):
    ''' use the dictionary, and the pan_genome_reference fasta file
    to create four different files for each cluster with the
    sequences of the genes in each category
    These will be used to postprocess the rare genes, and later
    to compare between different clusters.'''
    print("Getting the gene sequences....")
    counts = {}
    outputs = {}
    types = ["core", "soft_core", "inter", "rare"]
    for type in types:
        outputs[type] = open(os.path.join(d, type + "_genes.fa"), "w")
        counts[type] = 0
    with open(os.path.join(d, "merged_pan_genome_reference.fa")) as handle:
        for values in SimpleFastaParser(handle):
            gene_name = values[0].split()[-1]
            protein_sequence = translate(values[1])
            if gene_freqs[gene_name] < 0.15:
                outputs["rare"].write(">" + values[0] + "\n" + protein_sequence + "\n")
                counts["rare"] += 1
            elif gene_freqs[gene_name] < 0.95:
                outputs["inter"].write(">" + values[0] + "\n" + protein_sequence + "\n")
                counts["inter"] += 1
            elif gene_freqs[gene_name] < 0.99:
                outputs["soft_core"].write(">" + values[0] + "\n" + protein_sequence + "\n")
                counts["soft_core"] += 1
            else:
                outputs["core"].write(">" + values[0] + "\n" + protein_sequence + "\n")
                counts["core"] += 1
    for o in outputs:
        outputs[o].close()

    with open(os.path.join(d, "merged_summary_statistics.txt"), "w") as out:
        for t in types:
            out.write(t + "\t" + str(counts[t]) + "\n")
        out.write("Total\t" + str(sum(counts.values())) + "\n")
    return



def run(args):
    # get the classification of genes to their cluster
    merge_files(args.d)
    classify_genes(args.d)
    return

def get_options():
    parser = argparse.ArgumentParser(description='Extract the gene sequences from roary outputs, and merge rare genes')
    # input options
    parser.add_argument('--d', required=False,
                        type=str, default =  "/lustre/scratch118/infgen/team216/gh11/e_coli_collections/poppunk/new_roary/9_1558536260",
                        help='path to inpqut directory [%(default)s]')
    return parser.parse_args()


if __name__ == "__main__":
    # get arguments from user
    options = get_options()
    run(options)
