import argparse
import os
from Bio.SeqIO.FastaIO import SimpleFastaParser
from csv import reader
import networkx as nx
import subprocess
import glob
import string
import random
import shutil
import gffutils
import time

def get_gene_reps(d, gff_jobs_file, makeblastdb, blastn, cpus, ident, l, fg, lg):
    ''' Get a list of all the gene reps for a gene cluster
     create an empty fasta file for all the reps'''
    print("Getting gene reps...")
    d = os.path.abspath(d)
    cluster = d.split("/")[-1].split("_")[0]
    out_presence_absence = open(os.path.join(d, "new" + str(fg) + "_gene_presence_absence.Rtab"), "w")
    out_ref = open(os.path.join(d, "new" + str(fg) + "_pan_genome_reference.fa"), "w")
    out_cp = open(os.path.join(d, "new" + str(fg) + "_clustered_proteins"), "w")
    cnt = 0

    member_to_gene = {} ## member -> gene
    curr_genomes = set() ## list of all the current genomes
    reps_per_genome = {}
    genes = set()
    tmp_out = "tmp_" + ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(10))
    if not os.path.exists(tmp_out):
        os.makedirs(tmp_out)

    with open(os.path.join(d, "gene_presence_absence.csv")) as f:
        for toks in reader(f):
            if toks[0] == "Gene":
                rep_indexes = toks.index("Avg group size nuc") + 1
                genomes = toks[rep_indexes:]
                out_presence_absence.write("Gene\t" + "\t".join(genomes) + "\n")
                continue
            if cnt < fg:
                cnt += 1
                continue
            if cnt > lg:
                break
            gene_name = toks[0].replace("/", "_")
            genes.add(gene_name)
            cnt += 1 ## counter to stop after X number of genes
            print(gene_name)
            curr_reps = toks[rep_indexes:]
            for i in range(0, len(genomes)):
                curr_genome = genomes[i]
                curr_genome_reps = curr_reps[i]
                if curr_genome_reps == "":
                    continue
                if curr_genome not in reps_per_genome:
                    reps_per_genome[curr_genome] = []
                for r in curr_genome_reps.split("\t"):
                    reps_per_genome[curr_genome].append(r)
                    member_to_gene[r] = gene_name
    sep_genes(cluster, tmp_out, genomes, gff_jobs_file, genes, reps_per_genome, member_to_gene, makeblastdb, blastn, cpus, ident, l,  out_presence_absence, out_ref, out_cp)
    out_ref.close()
    out_presence_absence.close()
    out_cp.close()
    shutil.rmtree(tmp_out)
    return


def get_gene_sequences(tmp_out, curr_genome, gff_file, fasta_file, member_to_gene, reps_per_genome, gene_outs):
    '''Convert a gff file with the appended FASTA to protein/all fasta file
    gff_file = input gff file
    out_files = output files to write'''
    r = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(10))
    db = gffutils.create_db(gff_file, dbfn= r + '.db', force=True)
    curr_members = reps_per_genome[curr_genome]
    for member in curr_members:
        seq = db[member].sequence(fasta_file)
        curr_out = gene_outs[member_to_gene[member]]
        curr_out.write(">"+ curr_genome +"|" + member + "\n" + seq + "\n")
    os.remove(r + '.db')
    return


def blast_gene_cluster(tmp_out, gene_file, makeblastdb, blastn, cpus, ident, l):
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

    tmp_blast = os.path.join(tmp_out, ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(10)) + "_gene_blast_results.tab")
    subprocess.call([makeblastdb, "-in", gene_file, "-dbtype", "nucl"], stderr=subprocess.STDOUT)
    subprocess.call([blastn, "-query", gene_file,
                    "-db", gene_file, "-out", tmp_blast, "-num_threads",
                    str(cpus), "-dust", "no", "-outfmt", "6 qseqid sseqid pident length qlen slen evalue bitscore", "-evalue", "0.01"], stderr=subprocess.STDOUT)

    with open(tmp_blast) as f:  # add the edges to the graph
        for line in f:
            toks = line.strip().split()
            identity = float(toks[2])
            coverage = float(toks[3]) /  max(float(toks[4]), float(toks[5]))
            if coverage > l and identity > ident:
                G.add_edge(toks[0],toks[1])
    cc = nx.connected_components(G)
    files_to_delete = glob.glob(gene_file + "*")
    for f in files_to_delete:
        os.remove(f)
    os.remove(tmp_blast)
    return cc, rep_to_seq

def rewrite_roary_outputs(name, genomes, cc, rep_to_seq, out_presence_absence, out_ref, out_cp):
    ''' rewrite the roary outputs that matter which are:
    1. gene_presence_absence.Rtab
    2. pan_genome_reference.fasta = chooses longest rep and with least number of Ns
    3. clustered_proteins
    (4. gene_presence_absence.csv? -> maybe eventually if I need it, much more complicated)
    '''
    print("Rewriting outputs")
    new_genes = {}
    cnt = -1
    cc = list(cc)
    for c in cc:
        cnt += 1
        curr_genomes = []
        ## write reps for new cluster in clustered proteins output
        curr_name = name + "_" + str(cnt)
        out_cp.write(curr_name + ":")
        best_rep = {"name":"", "length":0, "Ns": 1000000, "seq":""}
        for rep in c:
        #    print("writing: " + rep)
            out_cp.write("\t" + rep)
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

def sep_genes(cluster, tmp_out, genomes, gff_jobs_file, genes, reps_per_genome, member_to_gene, makeblastdb, blastn, cpus, ident, l, out_presence_absence, out_ref, out_cp):
    ''' get the reps from all the GFF files and add them
    to the gene fasta for a single gene '''
    gene_outs = {}
    gene_file_names = {}
    for gene in genes:
        ## open a file with dictionary for out
        gene_file_names[gene] = os.path.join(tmp_out,''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(10)) + "_gene.fasta")
        gene_outs[gene] = open(gene_file_names[gene], "w")
    cnt = 0
    with open(gff_jobs_file) as f:
        for line in f:
            toks = line.strip().split("\t")
            if line.startswith("ID"):
                gff_loc = toks.index("New_annot_loc")
                poppunk_loc = toks.index("Poppunk_cluster")
                assembly_loc = toks.index("Assembly_Location")
                continue
            if toks[poppunk_loc] != cluster:
                continue
            gff_file = toks[gff_loc]
            curr_genome = gff_file.split("/")[-1].replace(".gff", "")
            if curr_genome not in reps_per_genome: ## this genome doesn't have this gene
                continue

            assemblies = toks[assembly_loc].split(",")
            chosen_assembly = None
            for a in assemblies:
                curr = a.split("/")[-1].split(".")[0]
                if curr.lower() in gff_file.lower():
                    chosen_assembly = a

            if chosen_assembly is None:
                print("gff: %s\nassemblies: %s" %(gff_file, str(assemblies)))

            cnt += 1
            print("Getting genes from genome: %d" %cnt)
            print("gff: %s\nassemblies: %s" %(gff_file, str(chosen_assembly)))
            start = time.time()
            get_gene_sequences(tmp_out, curr_genome, gff_file, chosen_assembly, member_to_gene, reps_per_genome, gene_outs)
            end = time.time()
            print(end - start)
    ## close all the current gene files
    for gene in gene_outs:
        gene_outs[gene].close()

    ## write the output
    for gene in gene_file_names:
        cc, rep_to_seq = blast_gene_cluster(tmp_out, gene_file_names[gene], makeblastdb, blastn, cpus, ident, l)
        rewrite_roary_outputs(gene, genomes, cc, rep_to_seq, out_presence_absence, out_ref, out_cp)
    return


def run(args):
    # get the classification of genes to their cluster
    get_gene_reps(args.d, args.metdata_file, args.makeblastdb, args.blastn, args.cpus, args.ident, args.l, args.fg, args.lg)
    return

def get_options():
    parser = argparse.ArgumentParser(description='Extract the gene sequences from roary outputs, and merge rare genes')
    # input options
    parser.add_argument('--d', required=False,
                        type=str, default =  "/lustre/scratch118/infgen/team216/gh11/e_coli_collections/poppunk/new_roary/48_1558538033",
                        help='path to inpqut directory [%(default)s]')
    parser.add_argument('--metdata_file', required=False,
                            type=str, default =  "/lustre/scratch118/infgen/team216/gh11/e_coli_collections/FILTERED_MD_FINAL_ALL.tab",
                            help='File with the list of GFF files and poppunk clusters[%(default)s]')
    parser.add_argument('--blastn',
                        type=str, default = "/software/pubseq/bin/ncbi_blast+/blastn",
                        help='blastn executable [%(default)s]')
    parser.add_argument('--makeblastdb',
                        type=str, default = "/software/pubseq/bin/ncbi_blast+/makeblastdb",
                        help='makeblastdb executable [%(default)s]')
    parser.add_argument('--cpus',
                        type=int, default = 4,
                        help='Number of CPUs to use [%(default)s]')
    parser.add_argument('--ident',
                        type=float, default = 80,
                        help='identity threshold to use for blastn [%(default)s]')
    parser.add_argument('--l',
                        type=float, default = 0.90,
                        help='length cutoff to use for sequence alignments [%(default)s]')
    parser.add_argument('--fg',
                        type=int, default = 0,
                        help='First gene to postprocess [%(default)s]')
    parser.add_argument('--lg',
                        type=int, default = 1,
                        help='Last gene to postprocess [%(default)s]')
    return parser.parse_args()


if __name__ == "__main__":
    # get arguments from user
    options = get_options()
    run(options)
