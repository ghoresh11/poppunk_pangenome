import argparse
import os
from csv import reader
import numpy as np
import operator
from Bio.Seq import reverse_complement


def get_gene_reps(d):
    ''' Get a list of all the gene reps for a gene cluster
     create an empty fasta file for all the reps'''
    print("Getting gene reps...")
    cluster = d.split("/")[-1].split("_")[0]
    reps = {}  ## genome -> rep gene -> gene cluster
    with open(os.path.join(d, "gene_presence_absence.csv")) as f:
        for toks in reader(f):
            if toks[0] == "Gene":
                rep_indexes = toks.index("Avg group size nuc") + 1
                genomes = toks[rep_indexes:]
                for g in genomes:
                    reps[g] = {}
                continue
            gene_name = toks[0].replace("/", "_")
            curr_reps = toks[rep_indexes:]
            for g in genomes:
                for r in curr_reps:
                    if r == "":
                        continue
                    reps[g][r] = gene_name
    return reps, genomes


def gff_to_fasta(gff_files, gff_file, reps, gene_lengths):
    '''Convert a gff file with the appended FASTA to protein/all fasta file
    gff_file = input gff file
    out_files = output files to write'''
    curr_genome = gff_file.split("/")[-1].replace(".gff", "")
    gff_files[curr_genome] = gff_file
    print(curr_genome)
    curr_reps = reps[curr_genome]
    contigs = {}
    with open(gff_file) as f:
        for line in f:
            if line.startswith("##FASTA"):
                break
            if line.startswith("#"):
                continue
            toks = line.strip().split("\t")
            if toks[2] != "CDS":
                continue
            name = toks[-1].split(";")[0].replace("ID=","")
            if name not in curr_reps: ## this CDS does not enter the roary ouput
                continue
            ## get the length of the gene
            gene_length = int(toks[4]) - int(toks[3]) + 1

            curr_gene = curr_reps[name]
            if curr_gene not in gene_lengths:
                gene_lengths[curr_gene] = {}
            gene_lengths[curr_gene][curr_genome + "|" + name] =  gene_length
    return

def read_gffs(gff_jobs_file, gene_reps):
    ''' get the reps from all the GFF files and add them
    to the gene fasta'''
    print("Getting gene sequences from GFF files...")
    gene_lengths = {} ## keep track of the length of each rep for each gene
    gff_files = {}
    with open(gff_jobs_file) as f:
        for line in f:
            line = line.strip()
            gff_to_fasta(gff_files, line, gene_reps, gene_lengths)
    return gff_files, gene_lengths


def sep_genes_by_size(gene_lengths, length_diff):
    ''' seperate the gene clusters according to their size
    So genes that vary by more than 0.2 in length are not in
    the same cluster (split them into two clusters)'''
    cc = {}
    for gene in gene_lengths:
        lengths_members = {}
        for member in gene_lengths[gene]:
            curr_length = gene_lengths[gene][member]
            min_curr_length = curr_length - length_diff * curr_length
            max_curr_length = curr_length + length_diff * curr_length
            chosen_key = None
            for l in lengths_members:
                if min_curr_length <= l and l <= max_curr_length: ## l is within the range:
                    lengths_members[l][member] = l
                    chosen_key = l
                    break
            if chosen_key is not None:
                ## update length to be a weighted mean of mean so far and new length
                num_prev_members = len(lengths_members[l].keys()) - 1
                new_length = int(np.mean([chosen_key] * num_prev_members + [curr_length]))
                if new_length == chosen_key:
                    continue
                lengths_members[new_length] = lengths_members[chosen_key]
                del lengths_members[chosen_key]
            else:
                lengths_members[curr_length] = {}
                lengths_members[curr_length][member] = curr_length
        cc[gene] = lengths_members.values()
    return cc

def get_rep_sequence(gene_id, gff_file):
    ''' get the rep sequences from the gff file '''
    is_contig = False
    fasta = False
    contig_seq = ""
    contig = "$$$"
    with open(gff_file) as f:
        for line in f:
            if is_contig and line.startswith(">"):
                break
            if is_contig:
                contig_seq += line.strip()
            if line.startswith(">" + contig):
                is_contig = True
                continue
            if fasta:
                continue
            if line.startswith("##FASTA"):
                fasta = True
                continue
            if line.startswith("#"):
                continue
            toks = line.strip().split("\t")
            if toks[2] != "CDS":
                continue
            name = toks[-1].split(";")[0].replace("ID=","")
            if name != gene_id: ## this CDS does not enter the roary ouput
                continue
            contig = toks[0]
            start = int(toks[3]) - 1
            stop = int(toks[4])
            strand = toks[6]
    seq = contig_seq[start:stop]
    if strand == "-":
        seq = reverse_complement(seq)
    return seq

def rewrite_for_gene(name, genomes, cc, gff_files, out_presence_absence, out_ref, out_cp):
    ''' rewrite the roary outputs that matter which are:
    1. gene_presence_absence.Rtab
    2. pan_genome_reference.fasta = chooses longest rep and with least number of Ns
    3. clustered_proteins
    (4. gene_presence_absence.csv? -> maybe eventually if I need it, much more complicated)
    '''
    new_genes = {}
    cnt = -1
    for c in cc:
        cnt += 1
        curr_genomes = []
        ## write reps for new cluster in clustered proteins output
        curr_name = name + "_" + str(cnt)
        out_cp.write(curr_name + ":")
        for rep in c:
            out_cp.write("\t" + rep.split("|")[1])
            curr_genomes.append(rep.split("|")[0]) ## save genomes with this gene
        out_cp.write("\n")
        ## update presence absence file
        out_presence_absence.write(curr_name)
        for g in genomes:
            if g in curr_genomes:
                out_presence_absence.write("\t1")
            else:
                out_presence_absence.write("\t0")
        out_presence_absence.write("\n")
        rep = max(c.iteritems(), key=operator.itemgetter(1))[0]
        out_ref.write(">" + rep + " " + curr_name  +  "\n" )
        out_ref.write(get_rep_sequence(rep.split("|")[1], gff_files[rep.split("|")[0]]) + "\n")
    return


def rewrite_roary_outputs(d, genomes, cc, gff_files):
    ''' rewrite the roary outputs after splitting genes
    by length'''
    out_presence_absence = open(os.path.join(d, "new_gene_presence_absence.Rtab"), "w")
    out_presence_absence.write("Gene\t" + "\t".join(genomes) + "\n")
    out_ref = open(os.path.join(d, "new_pan_genome_reference.fa"), "w")
    out_cp = open(os.path.join(d, "new_clustered_proteins"), "w")
    for gene in cc:
        rewrite_for_gene(gene, genomes, cc[gene], gff_files, out_presence_absence, out_ref, out_cp)
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
    with open(os.path.join(d, "new_gene_presence_absence.Rtab")) as f:
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
    with open(os.path.join(d, "new_pan_genome_reference.fa")) as handle:
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

    with open(os.path.join(d, "new_summary_statistics.txt"), "w") as out:
        for t in types:
            out.write(t + "\t" + str(counts[t]) + "\n")
        out.write("Total\t" + str(sum(counts.values())) + "\n")
    return



def run(args):
    # get the classification of genes to their cluster
    gene_reps, genomes = get_gene_reps(args.d)
    gff_files, gene_lengths = read_gffs(args.g, gene_reps)
    cc = sep_genes_by_size(gene_lengths, args.l)
    rewrite_roary_outputs(args.d, genomes, cc, gff_files)
    classify_genes(args.d)
    return

def get_options():
    parser = argparse.ArgumentParser(description='Extract the gene sequences from roary outputs, and merge rare genes')
    # input options
    parser.add_argument('--d', required=False,
                        type=str, default =  "/lustre/scratch118/infgen/team216/gh11/e_coli_collections/poppunk/new_roary/51_1558533994",
                        help='path to inpqut directory [%(default)s]')
    parser.add_argument('--g', required=False,
                            type=str, default =  "/lustre/scratch118/infgen/team216/gh11/e_coli_collections/poppunk/new_roary/new_jobs_corrected/jobs_51.txt",
                            help='File with the list of GFF files [%(default)s]')
    parser.add_argument('--l',
                        type=float, default = 0.2,
                        help='length cutoff to use for sequence alignments [%(default)s]')
    return parser.parse_args()


if __name__ == "__main__":
    # get arguments from user
    options = get_options()
    run(options)
