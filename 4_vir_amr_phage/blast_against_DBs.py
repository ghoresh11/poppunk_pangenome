import argparse
import os
import subprocess


def build_blast_dbs(args):
    ''' create a blast DB for the DB files
    return: a dictionary with db_name -> db_fasta_file_loc'''
    print("Making the BLAST databases....")
    db_dir = "/lustre/scratch118/infgen/team216/gh11/DBs"
    dbs = {"phasta_bact" : "bacteria_all_select_060319.db",
    "virulence": "clean_pathotype_genes.db",
    "phasta_prophage": "prophage_virus_060319.db",
    "resfinder" : "resfinder_cleaned_060319.db",
    "card": "card_ariba_cleaned_080319.fa",
    "arg-annot": "srst_ARGannot_r2_080319_cleaned.fasta"}
    #dbs = {"resfinder" : "resfinder_cleaned_060319.db"}
    for db in dbs:
        makeblastdb_command = map(
        str, [args.makeblastdb, "-in", os.path.join(db_dir, dbs[db]), "-dbtype", "prot"])
        #subprocess.call(makeblastdb_command)
        dbs[db] = os.path.join(db_dir, dbs[db])
    return dbs

def get_input_dirs(input_dir, cluster):
    ''' check all directories in the input dir
    return a list of directories that have all the required files'''
    print("Getting the input directories...")
    input_dir = os.path.abspath(input_dir)
    directories = [x[0] for x in os.walk(input_dir)]
    dirs_to_return = []
    for d in directories:
        if os.path.isfile(os.path.join(d, "pan_genome_reference.fa")) and os.path.isfile(os.path.join(d, "gene_presence_absence.csv")) and os.path.isfile(os.path.join(d, "gene_presence_absence.Rtab")):
        #    for debugging, uncomment:
            if d.split("/")[-1].split("_")[0] != cluster:
                continue
            dirs_to_return.append(d)
    return dirs_to_return


def blast_all_clusters(input_dirs, dbs, args):
    ''' iterate over the directories and blast all the
    rare genes, core genes, ... against each DB'''
    print("BLASTING against all DBs...")
    for d in input_dirs:
        for gene_class in ["rare", "inter", "soft_core", "core"]:
            for db in dbs:
                print("Blast CLuster %s, Gene Type: %s, Against: %s" %(d, gene_class, db))
                query_file = os.path.join(d, gene_class + "_genes.fa")
                out_file = os.path.join(d, gene_class + "_" + db + "_blast_result.tab")
                blastp_command = map(str, [args.blastp, "-db", dbs[db],
                               "-query", query_file, "-out",
                               out_file,
                               "-outfmt", "6 qseqid sseqid pident length qlen slen evalue bitscore",
                                "-evalue", args.min_blast_evalue,
                                "-num_threads", args.c])
                subprocess.call(blastp_command)
    return

def count_num_genes(dirs):
    ''' get the number of each gene class for each
    cluster to add to the summary'''
    print("Counting number of genes in each category...")
    num_genes = {}
    for d in dirs:
        curr_cluster = d.split("/")[-1].split("_")[0]
        num_genes[curr_cluster] = {}
        with open(os.path.join(d, "summary_statistics.txt")) as f:
            for l in f:
                if l.startswith("Core"):
                    num_genes[curr_cluster]["core"] = int(l.strip().split()[-1])
                    continue
                if l.startswith("Soft"):
                    num_genes[curr_cluster]["soft_core"] = int(l.strip().split()[-1])
                    continue
                if l.startswith("Shell"):
                    num_genes[curr_cluster]["inter"] = int(l.strip().split()[-1])
                    continue
                if l.startswith("Cloud"):
                    num_genes[curr_cluster]["rare"] = int(l.strip().split()[-1])
    return num_genes

def parse_blast_input(input_file, args):
    ''' count exactly how many hits there are from this category'''
    hits = set()
    with open(input_file) as f:
        for line in f:
            toks = line.strip().split("\t")
            identity = float(toks[2])
            coverage = float(toks[3]) / float(toks[5])

            if identity < args.min_identity or coverage < args.length_coverage:
                continue

            ## don't count anything twice, therefore it's a set
            hits.add(toks[0])
    return len(hits)


def summarise_blast_outputs(dirs, dbs, num_genes, args):
    ''' create a summary file
    that looks like:
    Cluster, Gene_class, DB, Count, Proportion
    1, rare, virulence, 10, 10/100
    1, rare, AMR, 20, 20/100
    etc
    '''
    print("Summarising the BLAST outputs...")
    out = open(args.cluster + "_gene_DB_hits.csv", "w")
    out.write("Cluster, Gene_class, DB, Count, Proportion\n")
    for d in dirs:
        curr_cluster = d.split("/")[-1].split("_")[0]
        for gene_class in ["core", "soft_core", "inter", "rare"]:
            for db in dbs:
                input_file = os.path.join(d, gene_class + "_" + db + "_blast_result.tab")
                cnt = parse_blast_input(input_file, args)
                out.write(",".join(map(str, [curr_cluster, gene_class, db, cnt, float(cnt) / num_genes[curr_cluster][gene_class]])) + "\n")
    out.close()
    return


def run(args):
    dbs = build_blast_dbs(args)
    dirs = get_input_dirs(args.input_dir, args.cluster)
    blast_all_clusters(dirs, dbs, args)
    num_genes = count_num_genes(dirs)
    summarise_blast_outputs(dirs, dbs, num_genes, args)
    return


def get_options():
    parser = argparse.ArgumentParser(description='Count how many of the DB genes are found in each gene class')
    # input options
    parser.add_argument('-mbe', '--min_blast_evalue', type=float,
                        help='Minimum BLAST evalue to use for an edge in the sequence similarity network [%(default)s]', metavar='FLOAT', default=0.01)
    parser.add_argument('-mi', '--min_identity', type=int,
                        help='Minimum BLAST identity to use for an edge in the sequence similarity network [%(default)s]', metavar='INT', default=75)
    parser.add_argument('-lc', '--length_coverage', type=float,
                        help='Minimum alignment coverage when comparing two sequences [%(default)s]', metavar='FLOAT', default=0.75)
    parser.add_argument('--makeblastdb', type=str,
                        help='makeblastdb executable [%(default)s]', default="/software/pubseq/bin/ncbi_blast+/makeblastdb", metavar='STR')
    parser.add_argument(
        '--blastp', type=str, help='blastn executable [%(default)s]', default="/software/pubseq/bin/ncbi_blast+/blastp", metavar='STR')
    parser.add_argument(
        '--c', help='cpus to use when running blast', default=8, type=int)
    parser.add_argument(
        '--cluster', help='Which cluster to run this run on (1-39)', default="39", type=str)
    parser.add_argument('--input_dir', required=False, type=str, default =  ".", help='path to input directory [%(default)s]')
    return parser.parse_args()


if __name__ == "__main__":
    # get arguments from user
    options = get_options()
    run(options)
