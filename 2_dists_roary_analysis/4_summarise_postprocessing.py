
import os

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
        if os.path.isfile(os.path.join(d, "new_summary_statistics.txt")):
            dirs_to_return[cluster] = os.path.join(input_dir, os.path.join(d, "new_summary_statistics.txt"))
    return dirs_to_return


def get_cluster_sizes():
    sizes = {}
    with open("cluster_sizes_updated.csv") as f:
        for line in f:
            if line.startswith("Cluster"):
                continue
            toks = line.strip().split(",")
            sizes[toks[0]] = toks[1]
    return sizes

files = get_input_dirs("/lustre/scratch118/infgen/team216/gh11/e_coli_collections/poppunk/new_roary/")
sizes = get_cluster_sizes()

out = open("new_summary_per_cluster.csv", "w")
out.write("cluster,variable,count,size\n")

for cluster in files:
    with open(files[cluster]) as f:
        for line in f:
            toks = line.strip().split("\t")
            if toks[0] == "Total":
                continue
            out.write(",".join([cluster,toks[0],toks[1],sizes[cluster]]) + "\n")

out.close()


out.close()
