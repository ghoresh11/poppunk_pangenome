import argparse
import os
import operator


def read_cogs_file(cogs_file="cog2cat60174_26-feb-2019.csv"):
    ''' read the cogs file into a dictionary
    with COG -> symbol, cat'''
    COGs_to_symbol_cat = {}
    with open(cogs_file) as f:
        for line in f:
            if line.startswith("COG ID"):
                continue
            toks = line.strip().split(",")
            cat = "-".join(toks[2:])
            cat = cat.replace("\"", "")
            COGs_to_symbol_cat[toks[0]] = {"symbol": toks[1], "cat": cat}
    return COGs_to_symbol_cat


def annot_to_symbol(desc, word_count):
    ''' use the description annotation to get
    a better COG symbol and category'''
    desc = desc.lower()
    words = desc.split()
    for w in words:
        if w not in word_count:
            word_count[w] = 0
        word_count[w] += 1

    if "phage" in desc or "conjuga" in desc or "transposase" in desc or "tail" in desc or "terminase" in desc or "antitermination" in desc or "plasmid" in desc or "transposon" in desc:
        return "X"
    if "hypothetical protein" in desc or "duf" in desc or desc == "protein":
        return "S"
    if "membrane" in desc:
        return "M"
    if "resistance" in desc or "multidrug" in desc or "crispr" in desc or "toxin" in desc or "effector" in desc or "immunity" in desc:
        return "V"
    if "secretion" in desc or "efflux" in desc:
        return "U"

    for i in range(1, 10):
        if "t" + str(i) + "ss" in desc:
            return "U"

    if "ribosom" in desc or "trna" in desc:
        return "J"
    if "rna polymerase" in desc or "transcriptional" in desc:
        return "K"
    if "protease" in desc:
        return "O"
    if "fimbria" in desc:
        return "N"
    if "dna polymerase " in desc or "replication" in desc:
        return "L"
    words = desc.split()
    for w in words:
        if (len(w) == 4 or w.endswith("-like")) and w.startswith("tra"):
            return "X"
        if w.endswith("ase"):
            return "R"
    if "transport" in desc:
        return "R"
    return "S"


def counts_to_summary(out, cluster_id, COG_counts, COGs_to_symbol_cat, gene_class, geneID_to_annot, word_count):
    '''create an output file with
    gene_name, COGs, Symbols, Categories
    This output can later be used in all downstream
    analysis for the functions of each cluster'''
    out = open(os.path.join(out, cluster_id + ".csv"), "w")
    out.write("Gene_name, COG, Symbol, Category, COG_Freq, Gene_Freq, Class\n")
    for gene_name in COG_counts:
        max_cog = max(COG_counts[gene_name].iteritems(),
                      key=operator.itemgetter(1))[0]
        size = sum(COG_counts[gene_name].values())
        freq = COG_counts[gene_name][max_cog] / float(size)
        symbol = "nd"
        if max_cog in COGs_to_symbol_cat:
            cat = COGs_to_symbol_cat[max_cog]["cat"]
            symbol = COGs_to_symbol_cat[max_cog]["symbol"]
        else:
            cat = geneID_to_annot[gene_name]
            symbol = annot_to_symbol(cat, word_count)
        out.write(",".join(map(str, [gene_name, max_cog, symbol, cat, freq,
                                     gene_class[gene_name]["freq"], gene_class[gene_name]["class"]])) + "\n")
    out.close()
    return


def summarise_COGs_per_cluster(geneID_to_COG, roary_dir):
    ''' return a dictionary with  cluster_name -> COG -> count'''
    COG_counts = {}
    clustered_protein_file = os.path.join(roary_dir, "clustered_proteins")
    with open(clustered_protein_file) as f:
        for line in f:
            toks = line.strip().split()
            gene_name = toks[0].replace(":", "")
            COG_counts[gene_name] = {}
            for member in toks[1:]:
                curr_cog = "nd"
                if member in geneID_to_COG:
                    curr_cog = geneID_to_COG[member]
                if curr_cog not in COG_counts[gene_name]:
                    COG_counts[gene_name][curr_cog] = 0
                COG_counts[gene_name][curr_cog] += 1
    return COG_counts


def gffs_from_jobs_file(jobs_file):
    ''' go over all the GFF files in a jobs file
    to get all the COGs'''
    geneID_to_COG = {}
    with open(jobs_file) as f:
        for line in f:
            gffs_to_COG(geneID_to_COG, line.strip())
    return geneID_to_COG


def gffs_to_COG(geneID_to_COG, gff_file):
    ''' add all lines of a GFF file to a
    dictionary with gene_id -> COG'''
    with open(gff_file) as f:
        for line in f:
            if line.startswith("##FASTA"):
                break
            if line.startswith("#"):
                continue
            toks = line.strip().split("\t")
            if toks[2] != "CDS":
                continue
            identifier = toks[-1]
            gene_name = identifier.split(";")[0].replace("ID=", "")
            COG = identifier.split("COG")
            if len(COG) == 1:  # if it's not in the dict it doesn't have a COG
                continue
            COG = COG[-1]
            COG = "COG" + COG.split(",")[0]
            COG = COG.split(";")[0]
            geneID_to_COG[gene_name] = COG
    return


def classify_genes(roary_dir):
    ''' go over the Rtab file for each roary cluster,
    calculate the frequency of each gene in the cluster
    return: dict classification of genes into core, soft_core, intermediate, rare'''
    gene_class = {}
    with open(os.path.join(roary_dir, "gene_presence_absence.Rtab")) as f:
        for line in f:
            if line.startswith("Gene"):
                continue
            toks = line.strip().split("\t")
            freq = sum(map(int, toks[1:])) / float(len(toks) - 1)
            if freq < 0.15:
                c = "rare"
            elif freq < 0.95:
                c = "inter"
            elif freq < 0.99:
                c = "soft_core"
            else:
                c = "core"
            gene_class[toks[0]] = {"freq": freq, "class": c}
    return gene_class


def summarise_functional_annotations(roary_dir):
    ''' Use the presence absence CSV file to count the functional annotations
    of the genes from the different groups in each cluster'''
    geneID_to_annot = {}
    with open(os.path.join(roary_dir, "gene_presence_absence.csv")) as f:
        for line in f:
            toks = line.strip().split(",")
            if toks[0] == "\"Gene\"":
                continue
            gene_name = toks[0].replace("\"", "")
            annotation = toks[2].replace("\"", "")
            geneID_to_annot[gene_name] = annotation
    return geneID_to_annot


def run(args):
    ''' read all the jobs files from the jobs dir
    return: a dictionary with all gene IDs to COGs'''
    COGs_to_symbol_cat = read_cogs_file()

    args.jobs_dir = os.path.abspath(args.jobs_dir)
    args.roary_dir = os.path.abspath(args.roary_dir)
    args.out = os.path.abspath(args.out)

    try:
        os.makedirs(args.out)
    except Exception:
        pass

    job_files = os.listdir(args.jobs_dir)
    roary_dirs = [x[0] for x in os.walk(args.roary_dir)]
    word_count = {}

    for f in job_files:
        if not f.startswith("jobs_"):
            continue
        # ncomment to test only on one small cluster
        if f != "jobs_9.txt":
            continue
        cluster_id = f.split("jobs_")[-1]
        cluster_id = cluster_id.replace(".txt", "")
        print("Getting functions for %s...." % cluster_id)
        geneID_to_COG = gffs_from_jobs_file(os.path.join(args.jobs_dir, f))
        for d in roary_dirs:
            if cluster_id + "_" in d:
                roary_dir = d
                break
        COG_counts = summarise_COGs_per_cluster(geneID_to_COG, roary_dir)
        geneID_to_annot = summarise_functional_annotations(roary_dir)
        gene_class = classify_genes(roary_dir)
        counts_to_summary(args.out, cluster_id, COG_counts,
                          COGs_to_symbol_cat, gene_class, geneID_to_annot, word_count)

    ## write out the word count -> to help with functional annotation
    with open(os.path.join(args.out, "word_count.csv"), "w") as out:
        for key in word_count:
            out.write(key + "," + str(word_count[key]) + "\n")
    return


def get_options():
    parser = argparse.ArgumentParser(
        description='Get the COG symbols and categories for all the genes from roary.')
    # input options
    parser.add_argument('--jobs_dir',
                        type=str, default="/lustre/scratch118/infgen/team216/gh11/e_coli_collections/poppunk/dists_analysis/gff_jobs/",
                        help='directory with all the GFF job files  [%(default)s]')
    parser.add_argument('--roary_dir',
                        required=False,
                        type=str, default="/lustre/scratch118/infgen/team216/gh11/e_coli_collections/poppunk/dists_analysis/roary_outputs/",
                        help='directory with the roary outs of all the job runs [%(default)s]')
    parser.add_argument('--out', help='Name of output directory [%(default)s]', default="functions/",
                        type=str,
                        required=False)
    return parser.parse_args()


if __name__ == "__main__":
    options = get_options()
    run(options)
    quit()
