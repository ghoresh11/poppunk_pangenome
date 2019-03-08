import argparse
import os
import subprocess
from csv import reader

def get_gff_locs(jobs_dir):
    ''' go over the gff_jobs directory and return
    a dictionary with:
    cluster -> genome identifier -> "gff_loc" -> full path to gff'''
    clusters = {}
    jobs_dir = os.path.abspath(jobs_dir)
    jobs_file = os.listdir(jobs_dir)
    for f in jobs_file:
        curr_cluster = f.split("_")[-1].replace(".txt", "")
        clusters[curr_cluster] = {}
        with open(os.path.join(jobs_dir, f)) as f_open:
            for line in f_open:
                basename = os.path.basename(line.strip())
                basename = basename.replace(".velvet.gff", "")
                basename = basename.replace(".gff", "")
                clusters[curr_cluster][basename] = {"gff_loc" : line.strip()}
    return clusters

def get_input_dirs(input_dir):
    ''' check all directories in the input dir
    return a list of directories that have all the required files'''
    print("Getting the input directories...")
    input_dir = os.path.abspath(input_dir)
    directories = [x[0] for x in os.walk(input_dir)]
    dirs_to_return = []
    for d in directories:
        if os.path.isfile(os.path.join(d, "pan_genome_reference.fa")) and os.path.isfile(os.path.join(d, "gene_presence_absence.csv")) and os.path.isfile(os.path.join(d, "gene_presence_absence.Rtab")):
            ## for debugging, uncomment:
            # if d.split("/")[-1].split("_")[0] != "39":
            #     continue
            dirs_to_return.append(d)
    return dirs_to_return


def get_gene_classes(input_dirs):
    ''' go over the Rtab file for each roary cluster,
    calculate the frequency of each gene in the cluster
    return: dict classification of genes into core, soft_core, intermediate, rare
    cluster -> gene -> class'''
    print("Calculating gene class....")
    gene_class = {}
    for d in input_dirs:
        curr_cluster = os.path.basename(d).split("_")[0]
        print(curr_cluster)
        gene_class[curr_cluster] = {}
        with open(os.path.join(d, "gene_presence_absence.Rtab")) as f:
            for line in f:
                if line.startswith("Gene"):
                    continue
                toks = line.strip().split("\t")
                freq = sum(map(int, toks[1:])) / float(len(toks)-1)
                if freq < 0.15:
                    gene_class[curr_cluster][toks[0]] = "rare"
                elif freq < 0.95:
                    gene_class[curr_cluster][toks[0]] = "inter"
                elif freq < 0.99:
                    gene_class[curr_cluster][toks[0]] = "soft_core"
                else:
                    gene_class[curr_cluster][toks[0]] = "core"
    return gene_class


def add_genes_to_genomes(clusters_genome_gene, cluster_gene_class, input_dirs):
    ''' go over the gene_presence_absence.csv file for each cluster.
    For each genome, add all its genes to the clusters_genome_gene dictionary
    and mark whether the gene is rare/core/inter/soft_core
    for rewriting the GFF files'''
    print("Getting class for each gene in genome...")
    for d in input_dirs:
        curr_cluster = os.path.basename(d).split("_")[0]
        print(curr_cluster)
        curr_genes = cluster_gene_class[curr_cluster]
        curr_genomes = clusters_genome_gene[curr_cluster]
        with open(os.path.join(d, "gene_presence_absence.csv")) as f:
            for toks in reader(f):
                if toks[0] == "Gene":
                    continue
                curr_gene = toks[0]
                members = toks[14:]
                for m in members:
                    if m == "":
                        continue
                    m = m.split("\t") ## a genome has multiple copies of this gene
                    for m2 in m:
                        genome_name = "_".join(m2.split("_")[:-1])
                        curr_genomes[genome_name][m2] = curr_genes[curr_gene]
    return


def rewrite_gff_files(input_dir, clusters_genome_gene_class):
    ''' for each genome, write a gff file with only CDSs of a specific class
    output four jobs files for running roary'''
    print("Rewriting GFF files...")
    input_dir = os.path.abspath(input_dir)
    out_dir = os.path.join(input_dir, "roary_on_class")
    gene_classes = ["rare", "inter", "soft_core", "core"]
    try:
        os.makedirs(out_dir)
        for gene_class in gene_classes:
            os.makedirs(os.path.join(out_dir, gene_class))
    except Exception:
        pass


    for cluster in clusters_genome_gene_class:
        print(cluster)
        for genome in clusters_genome_gene_class[cluster]:
            # if cluster != "39":
            #     continue

            curr_genome = clusters_genome_gene_class[cluster][genome]
            curr_gff_file = curr_genome["gff_loc"]
            gff_basename = os.path.basename(curr_gff_file)
            outputs = {}
            for gene_class in gene_classes:
                outputs[gene_class] = open(os.path.join(out_dir, gene_class, gff_basename), "w")

            with open(curr_gff_file) as f:
                    fasta = False
                    for line in f:
                        if line.startswith("##FASTA"):
                            fasta = True
                        if fasta:
                            for o in outputs:
                                outputs[o].write(line)
                            continue
                        if line.startswith("#"):
                            continue
                        toks = line.strip().split("\t")
                        if toks[2] != "CDS":
                            continue
                        gene_name = toks[-1].split(";")[0].replace("ID=","")
                        if gene_name not in curr_genome: ## this gene isn't in the csv file
                            ## for some reason a large porportion of genes
                            ## don't show up in the gene presence absence file at all
                            continue
                        gene_class = curr_genome[gene_name]
                        outputs[gene_class].write(line)

            for o in outputs:
                outputs[o].close()
    return (out_dir, gene_classes)


def run_roary(args, out_dir, gene_classes):
    ''' run roary on each gene class '''
    print("Running Roary")
    for gene_class in gene_classes:
        curr_dir = os.path.join(out_dir, gene_class)
        gff_files = os.listdir(curr_dir)
        gff_files = [os.path.join(curr_dir, f) for f in gff_files]
        job_name = gene_class + "_roary"
        mem = "10000"
        threads = "16s"

        lsf_prefix = ["bsub", "-q", "parallel", "-J", job_name, "-G", "team216","-o", job_name + ".o",
         "-e", job_name + ".e", '-R"select[mem>' + mem + '] rusage[mem='+ mem + '] span[hosts=1]"', '-M' + mem, "-n" + threads]

        command = map(str,["roary", "-p", threads, "-f", curr_dir, "-e", "-mafft", "-b", args.blastp, "-m", args.makeblastdb, "-s"]+ gff_files)
        ## submit the job
        subprocess.call(lsf_prefix + command)

        ## job_name=rare
        ## bsub -q parallel -J ${job_name} -G team216 -o ${job_name}.o -e ${job_name}.e -R"select[mem>10000]  rusage[mem=10000] span[hosts=1]" -M 10000 -n 16 roary -p 16 -e -mafft -b /software/pubseq/bin/ncbi_blast+/blastp -m /software/pubseq/bin/ncbi_blast+/makeblastdb -s *.gff
    return


def run(args):
    clusters_genome_gene = get_gff_locs(args.jobs_dir)
    input_dirs = get_input_dirs(args.input_dir)
    cluster_gene_class = get_gene_classes(input_dirs)
    add_genes_to_genomes(clusters_genome_gene, cluster_gene_class, input_dirs)
    out_dir, gene_classes = rewrite_gff_files(args.input_dir, clusters_genome_gene)
    run_roary(args, out_dir, gene_classes)
    return


def get_options():
    parser = argparse.ArgumentParser(
        description='Extract genes from each class in a GFF file and run roary on each class')
    # input options
    parser.add_argument('--jobs_dir', required=False, type=str,
                        default="/lustre/scratch118/infgen/team216/gh11/e_coli_collections/poppunk/dists_analysis/gff_jobs",
                        help='path to gff jobs directory [%(default)s]')
    parser.add_argument('--input_dir', required=False,
                        type=str, default =  ".",
                        help='path to input directory [%(default)s]')
    parser.add_argument('--blastp',
                        type=str, default = "/software/pubseq/bin/ncbi_blast+/blastp",
                        help='blastp executable [%(default)s]')
    parser.add_argument('--makeblastdb',
                        type=str, default = "/software/pubseq/bin/ncbi_blast+/makeblastdb",
                        help='makeblastdb executable [%(default)s]')
    return parser.parse_args()


if __name__ == "__main__":
    # get arguments from user
    options = get_options()
    run(options)
