import os

''' count how many core, soft-core, inter and rare genes there are in each cluster
and compare those numbers between roary and panaroo'''

def get_cluster_sizes():
    sizes = {}
    with open("/lustre/scratch118/infgen/team216/gh11/e_coli_collections/poppunk/cluster_sizes_updated.csv") as f:
        for line in f:
            if line.startswith("Cluster"):
                continue
            toks = line.strip().split(",")
            sizes[toks[0]] = toks[1]
    return sizes

def get_input_dirs(input_dir):
    ''' check all directories in the input dir
    return a list of directories that have all the required files'''
    input_dir = os.path.abspath(input_dir)
    directories = [d for d in os.listdir(input_dir)]
    dirs_to_return = []
    for d in directories:
        d = os.path.join(input_dir, d)
        cluster = d.split("/")[-1].split("_")[0]
        if not os.path.isdir(d):
            continue
        if os.path.isfile(os.path.join(d, "summary_statistics.txt")):
            dirs_to_return.append(os.path.join(input_dir, d))
    return dirs_to_return


print("cluster,tool,gene_type,count,size")
cluster_sizes = get_cluster_sizes()
roary_dirs = get_input_dirs("/lustre/scratch118/infgen/team216/gh11/e_coli_collections/poppunk/new_roary")
for d in roary_dirs:
    cluster = d.split("/")[-1].split("_")[0]
    with open(os.path.join(d, "summary_statistics.txt")) as f:
        for line in f:
            toks = line.strip().split("\t")
            if toks[0] == "Core genes":
                gene_type = "core"
            elif toks[0] == "Soft core genes":
                gene_type = "soft_core"
            elif toks[0] == "Shell genes":
                gene_type = "inter"
            elif toks[0] == "Cloud genes":
                gene_type = "rare"
            else:
                continue
            print(cluster + ",roary," + gene_type + "," + toks[-1] + "," + cluster_sizes[cluster])
    if not os.path.isfile(os.path.join(cluster, "gene_presence_absence.Rtab")):
        continue
    curr_file = os.path.join(cluster, "gene_presence_absence.Rtab")
    curr_counts = {"core": 0, "soft_core": 0, "inter": 0, "rare": 0}
    with open(curr_file) as f:
        for line in f:
            if line.startswith("Gene"):
                continue
            toks = line.strip().split()[1:]
            freq = sum(map(int, toks)) / float(len(toks))
            if freq > 0.99:
                curr_counts["core"] += 1
            elif freq > 0.95:
                curr_counts["soft_core"] += 1
            elif freq > 0.15:
                curr_counts["inter"] += 1
            else:
                curr_counts["rare"] += 1
    for o in curr_counts:
        print(cluster + ",panaroo," + o + "," + str(curr_counts[o]) + "," + cluster_sizes[cluster])
