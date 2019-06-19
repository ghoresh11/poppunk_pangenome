
import os
import subprocess


def get_input_dirs(input_dir):
    ''' check all directories in the input dir
    return a list of directories that have all the required files'''
    print("Getting the input directories...")
    input_dir = os.path.abspath(input_dir)
    directories = [d for d in os.listdir(input_dir)]
    dirs_to_return = {}
    for d in directories:
        d = os.path.join(input_dir, d)
        cluster = d.split("/")[-1].split("_")[0]
        if not os.path.isdir(d):
            continue
        if cluster in ["51"]: ## my test rabbits
            continue
        if os.path.isfile(os.path.join(d, "pan_genome_reference.fa")) and os.path.isfile(os.path.join(d, "gene_presence_absence.csv")) and os.path.isfile(os.path.join(d, "gene_presence_absence.Rtab")):
            dirs_to_return[cluster] = os.path.join(input_dir, d)
    return dirs_to_return


def get_cluster_sizes():
    sizes = {}
    with open("cluster_sizes_updated.csv") as f:
        for line in f:
            if line.startswith("Cluster"):
                continue
            toks = line.strip().split(",")
            sizes[toks[0]] = int(toks[1])
    return sizes

dirs = get_input_dirs("/lustre/scratch118/infgen/team216/gh11/e_coli_collections/poppunk/new_roary/")
sizes = get_cluster_sizes()

failed = [1]

for d in dirs:
    queue = "normal"
    cluster = d.split("/")[-1].split("_")[0]

    if int(cluster) not in failed:
        continue

    mem = str(sizes[cluster] * 10) ## the number of genomes determines the size I need of memory
                                    ## because i never run everything on more than 1000 genes at the same time
    job_name = "postprocess_" + cluster

    lsf_prefix = ["bsub", "-q", queue, "-J", job_name, "-G", "team216","-o", job_name + ".o",
         "-e", job_name + ".e", '-R"select[mem>' + mem + '] rusage[mem='+ mem + ']"', '-M' + mem]

    command = map(str,lsf_prefix + ["python", "4_post_process_roary.py",
    "--d", dirs[d], "--g",
    "/lustre/scratch118/infgen/team216/gh11/e_coli_collections/poppunk/new_roary/new_jobs_corrected/jobs_" + cluster + ".txt"])
    subprocess.call(command)
