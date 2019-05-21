import os
import subprocess


def get_input_dirs(input_dir):
    ''' check all directories in the input dir
    return a list of directories that have all the required files'''
    print("Getting the input directories...")
    input_dir = os.path.abspath(input_dir)
    directories = [x[0] for x in os.walk(input_dir)]
    dirs_to_return = []
    for d in directories:
        if os.path.isfile(os.path.join(d, "pan_genome_reference.fa")) and os.path.isfile(os.path.join(d, "gene_presence_absence.csv")) and os.path.isfile(os.path.join(d, "gene_presence_absence.Rtab")):
            # for debugging, uncomment:
            # if d.split("/")[-1].split("_")[0] != "39":
            #     continue
            dirs_to_return.append(d)
    return dirs_to_return


gene_types = ["rare", "soft_core"]
input_dirs = get_input_dirs("dists_analysis/roary_outputs/")

for d in input_dirs:
    for gene in gene_types:
        fasta_file = os.path.join(d, gene + "_genes.fa")
        cluster = d.split("/")[-1].split("_")[0]

        if cluster != "21":
            continue
        job_name = cluster + "_" + gene
        queue = "normal"
        mem = "20000"
        threads = "4"

        lsf_prefix = ["bsub", "-q", queue, "-J", job_name, "-G", "team216", "-o", job_name + ".o",
                      "-e", job_name + ".e", '-R"select[mem>' + mem + '] rusage[mem=' + mem + '] span[hosts=1]"', '-M' + mem, "-n" + threads]

        command = ["emapper.py", "-d", "bact", "--data_dir", "/lustre/scratch118/infgen/pathogen/pathpipe/eggnogmapper",
                   "--output", os.path.join("eggnog", job_name), "-m", "diamond", "-i", fasta_file, "--cpu", threads]

        #print(" ".join(lsf_prefix + command))
        subprocess.call(lsf_prefix + command)
