import os

## parse the roary outputs
roary_dir = "/lustre/scratch118/infgen/team216/gh11/e_coli_collections/poppunk/dists_analysis/roary_outputs/"
output_dirs = [x[0] for x in os.walk(roary_dir)]
## iterate over directories and see if a file exists
out = open("roary_summary_file.csv", "w")
out.write("cluster, core_genes, soft_core_genes, inter_genes, rare_genes, total_genes\n")
for od in output_dirs:
    input_file = os.path.join(od, "summary_statistics.txt")
    if os.path.isfile(input_file):
        cluster_id = od.split("/")[-1].split("_")[0]
        out_line = [cluster_id]
        with open(input_file) as f:
            for line in f:
                toks = line.strip().split()
                out_line.append(toks[-1])
        out.write(",".join(out_line) + "\n")
out.close()
