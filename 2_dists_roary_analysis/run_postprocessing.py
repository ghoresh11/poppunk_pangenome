
import os
import subprocess


MAX_GENES = 500

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
        # if cluster in ["51"]: ## my test rabbits
        #     continue
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

def get_num_genes(dirs):
    num_genes = {}
    for d in dirs:
        cnt = -1
        with open(os.path.join(dirs[d],"gene_presence_absence.Rtab")) as f:
            for line in f:
                cnt += 1
        num_genes[d] = cnt
    return num_genes


dirs = get_input_dirs("/lustre/scratch118/infgen/team216/gh11/e_coli_collections/poppunk/new_roary/")
sizes = get_cluster_sizes()
num_genes = get_num_genes(dirs)


failed = range(1,10)
# failed = [9]
## split the pipeline so it works on 5000 genes at a time?
## failed = range(11,51)
for d in dirs:
    queue = "normal"
    cluster = d.split("/")[-1].split("_")[0]

    if int(cluster) not in failed: ##change to 'in'
        continue

    mem = str(sizes[cluster] * 20) ## the number of genomes determines the size I need of memory
                                    ## because i never run everything on more than 1000 genes at the same time
    for i in range(0,num_genes[d],MAX_GENES):
        first = i
        last = i + MAX_GENES - 1
        if last > num_genes:
            last = num_genes

        job_name = "postprocess_" + str(i) + "_" + cluster

        lsf_prefix = ["bsub", "-q", queue, "-J", job_name, "-G", "team216","-o", job_name + ".o",
             "-e", job_name + ".e", '-R"select[mem>' + mem + '] rusage[mem='+ mem + ']"', '-M' + mem]

        command = map(str,lsf_prefix + ["python", "4_post_process_roary.py",
        "--d", dirs[d], "--g",
        "/lustre/scratch118/infgen/team216/gh11/e_coli_collections/poppunk/new_roary/new_jobs_corrected/jobs_" + cluster + ".txt",
        "--fg", str(first), "--lg", str(last)])
        subprocess.call(command)
