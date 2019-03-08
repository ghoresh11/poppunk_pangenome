import argparse
import os
import subprocess

def assemblies_to_reads(metadata_file):
    ''' parse the input metadata file so that each assembly
    has its equivalent read files.
    The READ files are required for running ARIBA.
    return: update dictionary of assembly -> Reads'''
    print("Getting READ filenames...")
    fa_to_reads = {}
    with open(metadata_file) as f:
        for line in f:
            toks = line.strip().split("\t")
            if toks[0] == "ID":
                assembly_index = toks.index("Assembly_Location")
                reads_index = toks.index("Reads_Location")
                continue
            assemblies = toks[assembly_index].split(",")
            reads = toks[reads_index].split(",")
            if len(reads) <= 2:
                for a in assemblies:
                    fa_to_reads[a] = reads
            else:
                for a in assemblies:
                    fa_to_reads[a] = []
                    basename = a.split("/")[-1]
                    basename = basename.split(".")[0]
                    for r in reads:
                        if basename in r:
                            fa_to_reads[a].append(r)
                    if a not in fa_to_reads:  # SANITY CHECK - SHOULDNT HAPPEN
                        print("NOT FOUND! \n assembly: %s\n gff: %s" %
                              (str(a), str(reads)))
    # for val in fa_to_reads:
    #     print("Key: %s, Values: %s" %(val, str(fa_to_reads[val])))
    return fa_to_reads


def get_clusters(clusters_file, fa_to_reads):
    ''' read the clusters_file so that for each assembly
    I know exactly which cluster its in
    return dict: assembly -> cluster'''
    clusters = {}
    with open(clusters_file) as f:
        for line in f:
            if line.startswith("Taxon"):
                continue
            toks = line.strip().split(",")
            if int(toks[1]) > 39:  # only looking at cluster 1-39
                continue
            name = os.path.basename(toks[0])
            clusters[name] = {
                "reads": fa_to_reads[toks[0]], "cluster": toks[1]}
    return clusters


def get_read_files(metadata_file, clusters):
    ''' add two read files to the assembly dict
    assembly -> cluster, read1, read2'''
    return


def run_ariba(clusters, out_dir):
    '''To run ariba:
    ## for every two files:
    ## It's probably best to create an ariba job for each cluster seperately on all
    ## genomes in that cluster
    ## Has to work with python3
    ariba run plasmidfinder_080319.prepareref reads1.fastq reads2.fastq out.run
    At the end, create a summary of all the output tsv files:
    ariba summary out.summary out.run1/report1.tsv out.run2/report2.tsv out.run3/report3.tsv --threads threads
    call a small ariba job 9000 times for all the genomes with very little memory and CPUs
    save the output file the name of the cluster + assembly so I run summary on all the
    ones that have the same prefix
    '''
    for item in clusters:
        # these have been long read sequenced, use a different method for them
        if len(clusters[item]["reads"]) != 2:
            continue

        out_path = os.path.join(
            out_dir, clusters[item]["cluster"] + "_" + item)

        job_name = clusters[item]["cluster"] + "_" + item
        mem = "500"
        threads = "1"
        lsf_prefix = ["bsub", "-J", job_name, "-G", "team216", "-o", job_name + ".o",
                      "-e", job_name + ".e", '-R"select[mem>' + mem + '] rusage[mem=' + mem + '] span[hosts=1]"', '-M' + mem, "-n" + threads]
        command = ["ariba", "run", "plasmidfinder_080319.prepareref", clusters[item]
                   ["reads"][0], clusters[item]["reads"][1], out_path, "--threads", threads]
        subprocess.call(lsf_prefix + command)
    return


def run(args):
    fa_to_reads = assemblies_to_reads(args.metadata_file)
    clusters = get_clusters(args.clusters_file, fa_to_reads)
    run_ariba(clusters, args.out_dir)
    return


def get_options():
    parser = argparse.ArgumentParser(
        description='Extract the gene sequences from roary outputs, and merge rare genes')
    # input options
    parser.add_argument('--metadata_file', required=False,
                        type=str, default="/lustre/scratch118/infgen/team216/gh11/e_coli_collections/FINAL_METADATA_CLEANED.csv",
                        help='path to metadata file [%(default)s]')
    parser.add_argument('--clusters_file', required=False,
                        type=str, default="/lustre/scratch118/infgen/team216/gh11/e_coli_collections/poppunk/poppunk_db_parallel/poppunk_db_parallel_clusters.csv",
                        help='path to clusters file [%(default)s]')
    parser.add_argument('--out_dir', required=False,
                        type=str, default="/lustre/scratch118/infgen/team216/gh11/e_coli_collections/plasmids/outputs/",
                        help='path to output directory[%(default)s]')

    return parser.parse_args()


if __name__ == "__main__":
    # get arguments from user
    options = get_options()
    run(options)
