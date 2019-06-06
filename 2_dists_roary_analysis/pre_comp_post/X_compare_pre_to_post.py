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
        if cluster in ["51"]: ## my test rabbits
            continue
        if os.path.isfile(os.path.join(d, "new_summary_statistics.txt")):
            dirs_to_return[cluster] = os.path.join(input_dir, d)
    return dirs_to_return


dirs = get_input_dirs(".")
out = open("pre_post_comparison.csv", "w")
out.write("cluster, pre_or_post, type, count\n")
for c in dirs:
    with open(os.path.join(dirs[c], "new_summary_statistics.txt")) as f:
        o = []
        for line in f:
            toks = line.strip().split("\t")
            out.write(c + ",post," + toks[0] + "," + toks[1] + "\n")
            o.append(toks[0])
    with open(os.path.join(dirs[c], "summary_statistics.txt")) as f:
        i = 0
        for line in f:
            toks = line.strip().split("\t")
            out.write(c + ",pre," + o[i] + "," + toks[2] + "\n")
            i += 1
out.close()
