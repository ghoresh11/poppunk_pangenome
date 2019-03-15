import os
import argparse
import csv

def sum_num_genes(gene_type):
    '''Count how many genes I start with'''
    cnt = 0
    num_rare_genes_in_cluster = {}
    with open("/Users/gh11/poppunk_pangenome/2_dists_roary_analysis/roary_summary_file.csv") as f:
        for line in f:
            if line.startswith("cluster"):
                continue
            toks = line.strip().split(",")
            if toks[1] != gene_type:
                continue
            cnt += int(toks[-1])
            if toks[0] not in num_rare_genes_in_cluster:
                num_rare_genes_in_cluster[toks[0]] = 0
            num_rare_genes_in_cluster[toks[0]] += int(toks[-1])
    return cnt, num_rare_genes_in_cluster


def check_orig_cluster():
    ''' read the assembly files and return a dictionary
    of assmbley -> cluster to be able to say later which gene cluster
    originated from which genome clusters'''
    genome_to_cluster = {}
    gff_job_files = os.listdir("gff_jobs/")
    for f in gff_job_files:
        curr_cluster = f.split("_")[-1].split(".")[0]
        with open(os.path.join("gff_jobs/", f)) as f_open:
            for line in f_open:
                curr_name = line.strip().split("/")[-1].replace(".velvet.gff","")
                curr_name = curr_name.replace(".gff","")
                genome_to_cluster[curr_name] = curr_cluster
    return genome_to_cluster




def read_clustered_proteins(gene_type, genome_to_cluster):
    ''' read the clustered proteins file into a dictionary to
    know what each protein cluster is built of'''
    pattern_count = {}
    num_members_count = {}
    with open(os.path.join(gene_type, "clustered_proteins")) as f:
        for line in f:
            toks = line.strip().split(":")
            name = toks[0]
            members = toks[1].strip().split("\t")
            cnt = {}
            for m in members:
                curr_genome = "_".join(m.split("_")[:-1])
                curr_cluster = genome_to_cluster[curr_genome]
                if curr_cluster not in cnt:
                    cnt[curr_cluster] = 0
                cnt[curr_cluster] += 1
            num_clusters = len(cnt)
            cluster_members = map(int,cnt.keys())
            cluster_members.sort()
            pattern = "-".join(map(str, cluster_members))

            if pattern not in pattern_count:
                pattern_count[pattern] = 0
            pattern_count[pattern] += 1

            if num_clusters not in num_members_count:
                num_members_count[num_clusters] = 0
            num_members_count[num_clusters] += 1
    return pattern_count, num_members_count


def generate_outputs(gene_type, pattern_count, num_members_count, sizes):
    ''' generate all relevant outputs to be loaded on to R'''
    matrix = [[0 for j in range(40)] for i in range(40)]

    matrix[0] = range(0,40)
    for i in range(0,40):
        matrix[i][0] = i

    with open(os.path.join(gene_type, "pattern_count.csv"), "w") as out:
        out.write("Pattern, Count, Num_members\n")
        for p in pattern_count:
            curr_members = map(int, p.split("-"))
            out.write(p + "," + str(pattern_count[p]) + "," + str(len(curr_members)) + "\n")
            for i in curr_members:
                for j in curr_members:
                    matrix[i][j] += 1

    # for i in range(1,40):
    #     for j in range(1,40):
    #         matrix[i][j] = matrix[i][j] / float(sizes[str(i)])

    with open(os.path.join(gene_type, "matrix.csv"), "wb") as out:
        writer = csv.writer(out)
        writer.writerows(matrix)

    with open(os.path.join(gene_type, "membership_count.csv"), "w") as out:
        out.write("Num_members, Count\n")
        for num_members in num_members_count:
            out.write(",".join(map(str, [num_members, num_members_count[num_members]])) + "\n")

    return

def run(args):
    num_orig_genes, num_rare_genes_in_cluster = sum_num_genes(args.type)
    genome_to_cluster = check_orig_cluster()
    pattern_count, num_members_count = read_clustered_proteins(args.type, genome_to_cluster)
    generate_outputs(args.type, pattern_count, num_members_count, num_rare_genes_in_cluster)
    return

def get_options():
    parser = argparse.ArgumentParser(description='Calculate the between and within cluster distances. Create files with roary jobs to run.')
    # input options
    parser.add_argument('--out', help='Name of output directory [%(default)s]', default = "dists_analysis",
                        type=str,
                        required=False)
    parser.add_argument('--type', help='What gene class is the analysis on [%(default)s]',
                        type=str,
                        default="rare")
   	# parser.add_argument('--merge_files', action='store_true', help='Merge the distances and clusters file into one', default=False)
    return parser.parse_args()


if __name__ == "__main__":
    # get arguments from user
    options = get_options()
    run(options)
