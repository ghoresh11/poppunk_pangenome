import argparse
import os
from numpy import mean, std
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.Seq import reverse_complement
from numpy import mean, std
import subprocess
import string
import random


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

def run_msa(curr_file):
    tmp_msa = curr_file + ".msa"
    out_msa = open(tmp_msa, "w")
    p = subprocess.Popen(["mafft", "--leavegappyregion", curr_file], stdout=out_msa, stderr=subprocess.PIPE)
    p.wait()
    out_msa.close()
    msa = {}
    with open(tmp_msa) as handle:
        for values in SimpleFastaParser(handle):
            msa[values[0]] = values[1]
            length = len(values[1])

    keys = list(msa.keys())
    mismatches = []
    for i in range(0, len(keys)-1):
        for j in range(i+1, len(keys)):
            if keys[i].split("|")[0] == keys[j].split("|")[0]:
                continue
            mismatches.append(sum(1 for x,y in zip(msa[keys[i]].lower(), msa[keys[j]].lower()) if x != y) / float(length))
    if len(mismatches) == 0:
        return None
    curr_mean = 1 - mean(mismatches)
    os.remove(tmp_msa)
    return curr_mean

def calc_variation(cluster1, cluster2, core_genes):
    ''' calculate the nucleotide variation within each cluster for each gene
    and the variation between every two clusters. Look at the difference between
    the two measures'''
    print("Calculating the pairwise variation...")
    ## within variation
    all_means = []
    for gene in core_genes:
        filenames = [os.path.join(str(cluster1), gene + ".fa"),os.path.join(str(cluster2), gene + ".fa")]
        tmp_both = "tmp_" + ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(5))
        with open(tmp_both, 'w') as outfile:
            for fname in filenames:
                with open(fname) as infile:
                    outfile.write(infile.read())
        curr_mean = run_msa(tmp_both)
        if curr_mean is None:
            continue
        all_means.append(curr_mean)
        os.remove(tmp_both)
    print("Mean: %f" %mean(all_means))
    print("SD: %f" %std(all_means))
    print("Min: %f" %min(all_means))
    return



def run(args):
    core_genes = get_core_genes()
    calc_variation(args.i, args.j, core_genes)
    return

def get_options():
    # need to run Roary with low cutoff to make this work on the reps I chose (60?)
    parser = argparse.ArgumentParser(
        description='Take all members of core genes, compare the variation within and between clusters to estimate identity cutoff.')
    parser.add_argument('-i', type=int, default = 1,
                            help='Cluster to run on [%(default)s]')
    parser.add_argument('-j', type=int, default = 2,
                            help='Cluster to run on [%(default)s]')
    return parser.parse_args()


if __name__ == "__main__":
    # get arguments from user
    options = get_options()
    run(options)
