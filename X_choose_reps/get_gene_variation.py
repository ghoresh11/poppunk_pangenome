import argparse
import os
from numpy import mean, std
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.Seq import reverse_complement
from numpy import mean, std
import subprocess

def get_core_genes():
    ''' read the gene presence absence file to get a list of all the core genes'''
    print("Creating a list of core genes....")
    core_genes = []
    with open("gene_presence_absence.Rtab") as f:
        for line in f:
            if line.startswith("Gene"):
                continue
            toks = line.strip().split("\t")
            freq = sum(map(int, toks[1:])) / float(len(toks)-1)
            if freq >= 0.95:
                core_genes.append(toks[0])
    return core_genes


def get_members_per_genome(core_genes):
    ''' read the clustered proteins file to get the members of the core genes
    in all the genomes'''
    print("Getting the members of all core genes...")
    gene_members = {}
    with open("clustered_proteins") as f:
        for line in f:
            toks = line.strip().split()
            gene_name = toks[0].replace(":","")
            if gene_name not in core_genes: ## not a core gene!
                continue
            for m in toks[1:]:
                gene_members[m] = gene_name
    return gene_members

def get_cluster_members():
    ''' assign each GFF file to a cluster '''
    print("Getting cluster and assembly info for all genomes...")
    genome_to_cluster = {}
    with open("300519_chosen.csv") as f:
        for line in f:
            toks = line.strip().split(',')
            if line.startswith("Name"):
                annot = toks.index("Annotation")
                assembly = toks.index("Assembly")
                cluster = toks.index("popppunk_cluster")
                continue
            annot_path = toks[annot].split("/")
            del annot_path[len(annot_path) - 2]
            annot_path = "/".join(annot_path)
            genome_to_cluster[toks[annot].split("/")[-1]] = {"cluster":toks[cluster], "assembly" : toks[assembly], "annot":annot_path}
    return genome_to_cluster

def create_files_per_cluster(i, core_genes, gene_members, genome_to_cluster):
    ''' for each gene, create a fasta file of all its members in each cluster
    cluster_gene.fasta'''
    outfiles = {}
    print("Generating gene fasta files for cluster %s..." %i)

    if not os.path.exists(i):
        os.makedirs(i)

    for gene in core_genes:
        outfiles[gene] = open(os.path.join(i, gene + ".fa"), "w")

    problem_files = set()
    gff_files = os.listdir(".")
    for f in gff_files:
        if not f.endswith(".gff"):
            continue
        curr_cluster = genome_to_cluster[f]["cluster"]
        if curr_cluster != i:
            continue
        curr_assembly = genome_to_cluster[f]["assembly"]
        contigs = {}
        with open(curr_assembly) as handle:
            for values in SimpleFastaParser(handle):
                contigs[values[0]] = values[1]

        with open(genome_to_cluster[f]["annot"]) as f_open:
            for line in f_open:
                if line.startswith("##FASTA"):
                    break
                if line.startswith("#"):
                    continue
                toks = line.strip().split("\t")
                name = toks[-1].split(";")[0].replace("ID=","")
                if name not in gene_members:
                    continue
                curr_gene = gene_members[name]
                curr_file = outfiles[curr_gene]
                if toks[0] not in contigs:
                    problem_files.add(f)
                    continue
                curr_seq = contigs[toks[0]][int(toks[3]) - 1 : int(toks[4])]
                if toks[6] == "-":
                    curr_seq = reverse_complement(curr_seq)
                curr_file.write(">" + str(curr_cluster) + "|" + name + "\n" + curr_seq + "\n")
    print(problem_files)
    for gene in outfiles:
        outfiles[gene].close()
    return


def run_msa(curr_file, between = False):
    out_msa = open("tmp_msa_file.fa", "w")
    p = subprocess.Popen(["mafft", "--leavegappyregion", curr_file], stdout=out_msa, stderr=subprocess.PIPE)
    p.wait()
    out_msa.close()
    msa = {}
    with open("tmp_msa_file.fa") as handle:
        for values in SimpleFastaParser(handle):
            msa[values[0]] = values[1]
            length = len(values[1])

    keys = list(msa.keys())
    mismatches = []
    for i in range(0, len(keys)-1):
        for j in range(i+1, len(keys)):
            if between and keys[i].split("|")[0] == keys[j].split("|")[0]:
                continue
            mismatches.append(sum(1 for x,y in zip(msa[keys[i]].lower(), msa[keys[j]].lower()) if x != y) / float(length))
    if len(mismatches) == 0:
        return None
    curr_mean = 1 - mean(mismatches)
    return curr_mean

def calc_variation(core_genes):
    ''' calculate the nucleotide variation within each cluster for each gene
    and the variation between every two clusters. Look at the difference between
    the two measures'''

    out = open("variation.csv", "w")
    out.write("Cluster1, Cluster2, Mean, SD\n")

    ## within variation
    for cluster in range(1,3):
        if cluster == 50:
            continue
        print("Calculating within cluster variation for cluster %d..." %cluster)
        all_means = []
        for gene in core_genes:
            curr_file = os.path.join(str(cluster), gene + ".fa")
            curr_mean = run_msa(curr_file)
            if curr_mean is None:
                continue
            all_means.append(curr_mean)
        out.write(",".join(map(str,[cluster, cluster, mean(all_means), std(all_means)])) + "\n")

    for cluster1 in range(1,2):
        for cluster2 in range(i+1,3):
            all_means = []
            for gene in core_genes:
                filenames = [os.path.join(str(cluster1), gene + ".fa"),os.path.join(str(cluster2), gene + ".fa")]
                with open('tmp_both', 'w') as outfile:
                    for fname in filenames:
                        with open(fname) as infile:
                            outfile.write(infile.read())
                curr_mean = run_msa(curr_file, between = True)
                if curr_mean is None:
                    continue
                all_means.append(curr_mean)
            out.write(",".join(map(str,[cluster1, cluster2, mean(all_means), std(all_means)])) + "\n")
    out.close()
    return


def run(args):
    core_genes = get_core_genes()
    gene_members = get_members_per_genome(core_genes)
    genome_to_cluster = get_cluster_members()
    for i in range(2,3):
        if i == 50:
            continue
        create_files_per_cluster(str(i), core_genes, gene_members, genome_to_cluster) # skip this if completed
    calc_variation(core_genes)
    return

def get_options():
    # need to run Roary with low cutoff to make this work on the reps I chose (60?)
    parser = argparse.ArgumentParser(
        description='Take all members of core genes, compare the variation within and between clusters to estimate identity cutoff.')
    return parser.parse_args()


if __name__ == "__main__":
    # get arguments from user
    options = get_options()
    run(options)
