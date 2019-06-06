import argparse
import os
from Bio.SeqIO.FastaIO import SimpleFastaParser
import subprocess
from Bio.Seq import reverse_complement
from Bio.Seq import translate
from csv import reader
import string
import random
import networkx as nx
import shutil
import glob

def get_gene_reps(d):
    ''' Get a list of all the gene reps for a gene cluster
     create an empty fasta file for all the reps'''
    print("Getting gene reps...")
    cluster = d.split("/")[-1].split("_")[0]
    tmp_out = cluster + "_tmp" + ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(5))
    os.makedirs(tmp_out)
    out_files = {}

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
            out_files[gene_name] = open(os.path.join(tmp_out, gene_name + ".fasta"), "w")
            curr_reps = toks[rep_indexes:]
            for g in genomes:
                for r in curr_reps:
                    if r == "":
                        continue
                    reps[g][r] = gene_name
    return reps, tmp_out, out_files, genomes


def gff_to_fasta(gff_file, reps, out_files, tmp_out):
    '''Convert a gff file with the appended FASTA to protein/all fasta file
    gff_file = input gff file
    out_files = output files to write'''
    curr_genome = gff_file.split("/")[-1].replace(".gff", "")
    print(curr_genome)
    curr_reps = reps[curr_genome]
    tmp_genome = os.path.join(tmp_out,"curr_genome.fa")
    out = open(tmp_genome, "w")
    contigs = {}
    with open(gff_file) as f:
        fasta = False
        for line in f:
            if fasta:
                out.write(line)
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
            if toks[0] not in contigs:
                contigs[toks[0]] = []
            contigs[toks[0]].append({"name": name, "start": int(toks[3]) - 1,
                                     "stop": int(toks[4]), "strand": toks[6]})
    out.close()
    # read the contigs and save the final fasta file
    with open(tmp_genome) as handle:
        for values in SimpleFastaParser(handle):
            curr_contig = values[0]
            if curr_contig not in contigs:  # no CDSs in this contig
                continue
            for cds in contigs[curr_contig]:
                if cds["name"] not in curr_reps:
                    continue
                out = out_files[curr_reps[cds["name"]]]
                out.write(">" + curr_genome +"|" + cds["name"] + "\n")
                seq = values[1][cds["start"]:cds["stop"]]
                if cds["strand"] == "-":
                    seq = reverse_complement(seq)
                out.write(seq + "\n")
    os.remove(tmp_genome)
    return

def read_gffs(gff_jobs_file, gene_reps, out_files, tmp_out):
    ''' get the reps from all the GFF files and add them
    to the gene fasta'''
    print("Getting gene sequences from GFF files...")
    with open(gff_jobs_file) as f:
        for line in f:
            line = line.strip()
            gff_to_fasta(line, gene_reps, out_files, tmp_out)
    for gene_file in out_files:
        out_files[gene_file].close()
        out_files[gene_file] = os.path.join(tmp_out, gene_file + ".fasta")
    return


def blast_gene_cluster(gene_file, tmp_out, makeblastdb, blastn, cpus, i, l):
    ''' Blast the rare genes of a file against each other
    return a dictionary with rare gene -> new rep'''
    G = nx.Graph()
    rep_to_seq = {}
    with open(gene_file) as handle:
        for values in SimpleFastaParser(handle):
            rep_to_seq[values[0]] =  values[1] ## this will provide num Ns and length
            G.add_node(values[0])
    if len(rep_to_seq.keys()) == 1: ## nothing to do
        os.remove(gene_file)
        return [rep_to_seq.keys()], rep_to_seq

    subprocess.call([makeblastdb, "-in", gene_file, "-dbtype", "nucl"], stderr=subprocess.STDOUT)
    subprocess.call([blastn, "-query", gene_file,
                    "-db", gene_file, "-out", os.path.join(tmp_out, "gene_blast_results.tab"), "-num_threads",
                    str(cpus), "-outfmt", "6 qseqid sseqid pident length qlen slen evalue bitscore", "-evalue", "0.01"], stderr=subprocess.STDOUT)
    with open(os.path.join(tmp_out, "gene_blast_results.tab")) as f:  # add the edges to the graph
        for line in f:
            toks = line.strip().split()
            identity = float(toks[2])
            coverage = float(toks[3]) /  max(float(toks[4]), float(toks[5]))
            if coverage > l and identity > i:
                G.add_edge(toks[0],toks[1])
    cc = nx.connected_components(G)
    files_to_delete = glob.glob(gene_file + "*")
    for f in files_to_delete:
        os.remove(f)
    return cc, rep_to_seq

def rewrite_roary_outputs(name, genomes, cc, rep_to_seq, out_presence_absence, out_ref, out_cp):
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
        best_rep = {"name":"", "length":0, "Ns": 1000000, "seq":""}
        for rep in c:
            out_cp.write("\t" + rep.split("|")[1])
            curr_genomes.append(rep.split("|")[0]) ## save genomes with this gene
            curr_seq = rep_to_seq[rep]
            curr_Ns = curr_seq.upper().count('N')
            ## choose the rep that has the least number of Ns and is the longest
            if len(curr_seq) > best_rep["length"] and  curr_Ns < best_rep["Ns"]:
                best_rep = {"name":rep, "length": len(curr_seq), "Ns": curr_Ns, "seq":curr_seq}
        out_cp.write("\n")
        ## keep consistency of chosen rep then name of gene
        out_ref.write(">" + best_rep["name"] + " " + curr_name  +  "\n" + best_rep["seq"] + "\n")
        ## update presence absence file
        out_presence_absence.write(curr_name)
        for g in genomes:
            if g in curr_genomes:
                out_presence_absence.write("\t1")
            else:
                out_presence_absence.write("\t0")
        out_presence_absence.write("\n")

    return


def blast_all_genes(d, out_files, genomes, tmp_out, makeblastdb, blastn, cpus, i, l):
    '''Blast all gene files
    Run blast on each gene file and split genes which now form more connected components'''
    out_presence_absence = open(os.path.join(d, "new_gene_presence_absence.Rtab"), "w")
    out_presence_absence.write("Gene\t" + "\t".join(genomes) + "\n")
    out_ref = open(os.path.join(d, "new_pan_genome_reference.fa"), "w")
    out_cp = open(os.path.join(d, "new_clustered_proteins"), "w")
    cnt = 0
    for gene in out_files:
        cc, rep_to_seq = blast_gene_cluster(out_files[gene], tmp_out, makeblastdb, blastn, cpus, i, l)
        rewrite_roary_outputs(gene, genomes, cc, rep_to_seq, out_presence_absence, out_ref, out_cp)
        # if cnt == 50:
        #     break
        cnt += 1
    out_ref.close()
    out_presence_absence.close()
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
    gene_reps, tmp_out, out_files, genomes = get_gene_reps(args.d)
    # read all the GFF files for this cluster, and save
    # all the gene reps to a temporary FASTA file
    read_gffs(args.gff_jobs_file, gene_reps, out_files, tmp_out)
    blast_all_genes(args.d, out_files, genomes, tmp_out, args.makeblastdb, args.blastn, args.cpu, args.i, args.l)
    shutil.rmtree(tmp_out) ## delete the temporary folder with all temp file
    # classify all the genes into core/rare etc and rewrite summary file
    classify_genes(args.d)
    return

def get_options():
    parser = argparse.ArgumentParser(description='Extract the gene sequences from roary outputs, and merge rare genes')
    # input options
    parser.add_argument('--d', required=False,
                        type=str, default =  "/lustre/scratch118/infgen/team216/gh11/e_coli_collections/poppunk/new_roary/51_1558533994",
                        help='path to inpqut directory [%(default)s]')
    parser.add_argument('--gff_jobs_file', required=False,
                            type=str, default =  "/lustre/scratch118/infgen/team216/gh11/e_coli_collections/poppunk/new_roary/new_jobs_corrected/jobs_51.txt",
                            help='File with the list of GFF files [%(default)s]')
    parser.add_argument('--blastn',
                        type=str, default = "/software/pubseq/bin/ncbi_blast+/blastn",
                        help='blastn executable [%(default)s]')
    parser.add_argument('--makeblastdb',
                        type=str, default = "/software/pubseq/bin/ncbi_blast+/makeblastdb",
                        help='makeblastdb executable [%(default)s]')
    parser.add_argument('--cpu',
                        type=int, default = 4,
                        help='Number of CPUs to use [%(default)s]')
    parser.add_argument('--i',
                        type=float, default = 95,
                        help='identity threshold to use for blastn [%(default)s]')
    parser.add_argument('--l',
                        type=float, default = 0.8,
                        help='length cutoff to use for sequence alignments [%(default)s]')
    return parser.parse_args()


if __name__ == "__main__":
    # get arguments from user
    options = get_options()
    run(options)
