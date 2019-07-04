
import os
import subprocess


MAX_GENES = 15

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


failed = [1]

for d in dirs:
    queue = "normal"
    cluster = d.split("/")[-1].split("_")[0]
    mem = "1000"
    if int(cluster) not in failed: ##change to 'in'
        continue

    for i in range(15, num_genes[d], MAX_GENES):
    #for i in range(0, failed[int(cluster)], MAX_GENES): ## only carry on until the number of genes that failed
        first = i
        last = i + MAX_GENES - 1
        if last > num_genes[d]:
            last = num_genes[d]

        job_name = "postprocess_" + cluster + "_" + str(i)
        print(job_name)
        lsf_prefix = ["bsub", "-q", queue, "-J", job_name, "-G", "team216","-o", job_name + ".o",
             "-e", job_name + ".e", '-R"select[mem>' + mem + '] rusage[mem='+ mem + ']"', '-M' + mem]

        command = map(str,lsf_prefix + ["python", "4_post_process_roary.py",
        "--d", dirs[d],
        "--fg", str(first), "--lg", str(last)])
        subprocess.call(command)
