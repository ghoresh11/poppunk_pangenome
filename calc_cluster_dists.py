import argparse
import os
from numpy import mean
import time

def assemblies_to_gffs(metadata_file):
    ''' parse the input metadata file so that each assembly
    has its equivalent GFF file.
    The GFF files are required for running ROARY.
    return: dictionary of assembly -> GFF'''
    print("Getting GFF filenames...")
    fa_to_gff = {}
    with open(metadata_file) as f:
        for line in f:
            toks = line.strip().split("\t")
            if toks[0] == "ID":
                assembly_index = toks.index("Assembly_Location")
                annot_index = toks.index("Annotation_Location")
                continue
            assemblies = toks[assembly_index].split(",")
            gffs = toks[annot_index].split(",")
            if len(gffs) == 1:
                for a in assemblies:
                    fa_to_gff[a] = gffs[0]
            else:
                for a in assemblies:
                    basename = a.split("/")[-1]
                    basename = basename.split(".")[0]
                    for g in gffs:
                        if basename in g:
                            fa_to_gff[a] = g
                    if a not in fa_to_gff: ## SANITY CHECK - SHOULDNT HAPPEN
                        print("NOT FOUND! \n assembly: %s\n gff: %s" % (str(a), str(gffs)))
    return fa_to_gff


def get_cluster_sizes(clusters_file, out):
    ''' go over the poppunk cluster output
    For each assembly, assign it a cluster membership
    Return 1) dictionary assembly -> cluster
        2) For each cluster, calculate its size. If it's smaller
    than min_cluster_size disregard it for ROARY by
    saving it in the cluster_size dictionary
    Output: a file with each clusters size '''
    print("Calculating cluster sizes...")
    try:
        os.makedirs(out)
    except Exception:
        pass
    clusters = {}
    cluster_size = {}
    with open(clusters_file) as f:
        for line in f:
            if line.startswith("Taxon"):
                continue
            toks = line.strip().split(",")
            if toks[1] not in cluster_size:
                cluster_size[toks[1]] = 0
            cluster_size[toks[1]] += 1
            clusters[toks[0]] = toks[1]
    ## generate cluster size output
    with open(os.path.join(out, "cluster_sizes.csv"),"w") as out:
        out.write("Cluster, Size\n")
        for cluster in cluster_size:
            out.write(cluster + "," + str(cluster_size[cluster]) + "\n")
    return (clusters, cluster_size)


def metadata_per_cluster(metadata_file, cluster_sizes, clusters, min_cluster_size, out):
    ''' generate a new clusters output so that I can generate bar/piecharts
    where I know exactly how many of each ST, pathotype, etc. are in each of the
    poppunk clusters
    output: melted dataframe with cluster and count of each number from each category
    return: nothing'''

    ## variables to ouput in file
    variables = ["Pathotype", "Country", "Continent", "Year","Isolation", "ST", "MASH"]
    clusters_metadata = {}
    with open(metadata_file) as f:
        for line in f:
            toks = line.strip().split("\t")
            if line.startswith("ID"):
                indexs = [toks.index(i) for i in variables]
                assembly_index = toks.index("Assembly_Location")
                continue
            assemblies = toks[assembly_index].split(",")
            for a in assemblies:
                if a not in clusters: ## some are problematic
                    continue
                cluster = clusters[a]
                if cluster_sizes[cluster] < min_cluster_size:
                    continue
                if cluster not in clusters_metadata:
                    clusters_metadata[cluster] = {}
                    for v in variables:
                        clusters_metadata[cluster][v] = {}
                for i in range(len(variables)):
                    curr_value = toks[indexs[i]]
                    v = variables[i]
                    ## this value has never been seen for this cluster
                    if curr_value not in clusters_metadata[cluster][v]:
                        clusters_metadata[cluster][v][curr_value] = 0
                    ## add count by one for this value in this variable in this cluster
                    clusters_metadata[cluster][v][curr_value] += 1

    ## generate the melted output file -> easy to work with in R
    out = open(os.path.join(out,"metadata_per_cluster.csv"),"w")
    out.write("cluster\tvariable\tvalue\tcount\n")
    for cluster in clusters_metadata:
        for var in clusters_metadata[cluster]:
            for value in clusters_metadata[cluster][var]:
                out.write("\t".join([cluster, var, value, str(clusters_metadata[cluster][var][value] / float(cluster_sizes[cluster]))]) + "\n")
    out.close()
    return

def get_dist_within_cluster(clusters, cluster_sizes, dists_file, out, min_cluster_size, fa_to_gff):
    ''' Due to memory constrains, read the dists file line by line.
    Check membership of both clusters, if they are different ignore,
    if they are the same but the cluster is too small also ignore
    (i.e. skip MOST lines of any calculation)
    Otherwise, save the core and acc dists of this cluster
    to be used for ROARY input and save all the GFF files of this clusters to
    a file to be used by ROARY
    Return: dictionary with cluster -> core and accessory distance'''
    print("Calculating within-cluster distances...")
    within_cluster_dist = {}
    gffs = {}
    cnt = 0
    with open(dists_file) as f:
        for line in f:
            toks = line.strip().split("\t")
            if line.startswith("Query"):
                continue
            if clusters[toks[0]] != clusters[toks[1]] or cluster_sizes[clusters[toks[0]]] < min_cluster_size:
                continue
            cluster = clusters[toks[0]]
            if cluster not in within_cluster_dist:
                within_cluster_dist[cluster] = {"core":[], "acc":[]}
                gffs[cluster] = set()
            within_cluster_dist[cluster]["core"].append(float(toks[2]))
            within_cluster_dist[cluster]["acc"].append(float(toks[3]))
            gffs[cluster].add(fa_to_gff[toks[0]])
            gffs[cluster].add(fa_to_gff[toks[1]])
            # if cnt == 10000: ## if testing uncomment here
            #     break
            cnt += 1

    ## calc the averages
    ## generate output so that this doesn't need to be repeated
    with open(os.path.join(out, "within_cluster_dist.csv"),"w") as f_out:
        f_out.write("Cluster, Core, Core_max, Acc, Acc_max\n")
        for cluster in within_cluster_dist:
            core_max = max(within_cluster_dist[cluster]["core"])
            acc_max = max(within_cluster_dist[cluster]["acc"])
            within_cluster_dist[cluster]["core"] = mean(within_cluster_dist[cluster]["core"])
            within_cluster_dist[cluster]["acc"] =  mean(within_cluster_dist[cluster]["acc"])
            f_out.write(",".join(map(str, [cluster, within_cluster_dist[cluster]["core"], core_max, within_cluster_dist[cluster]["acc"], acc_max])) + "\n")
    ## create the GFF outputs
    out = os.path.abspath(os.path.join(out,"gff_jobs"))
    try:
        os.makedirs(out)
    except Exception:
        pass
    for cluster in gffs:
        with open(os.path.join(out, "jobs_" + cluster + ".txt"), "w") as f_out:
            for gff in gffs[cluster]:
                f_out.write(gff + "\n")
    return within_cluster_dist


def add_clusters_to_dists(between_cluster_dist, cluster1, cluster2, core, accessory):
    ''' add the distance between two clusters to the between clusters dictionary '''
    if cluster1 not in between_cluster_dist:
        between_cluster_dist[cluster1] = {}
    if cluster2 not in between_cluster_dist[cluster1]:
        between_cluster_dist[cluster1][cluster2] = {"core":[], "accessory":[]}
    between_cluster_dist[cluster1][cluster2]["core"].append(core)
    between_cluster_dist[cluster1][cluster2]["accessory"].append(accessory)
    return

def between_cluster_dist(clusters, within_cluster_dist, dists_file, fa_to_gff, out):
    ''' read the dists file again, this time I only care about the
    distance between two clusters that I decided to keep, therefore
    the number of calculations on the whole file is still small '''
    print("Calculating between-cluster distances...")
    between_cluster_dist = {}
    with open(dists_file) as f:
        for line in f:
            toks = line.strip().split("\t")
            if line.startswith("Query"):
                continue
            if  clusters[toks[0]] == clusters[toks[1]] or clusters[toks[0]] not in within_cluster_dist or clusters[toks[1]] not in within_cluster_dist:
                continue
            cluster1 = int(clusters[toks[0]])
            cluster2 = int(clusters[toks[1]])
            if cluster1 < cluster2:
                add_clusters_to_dists(between_cluster_dist, cluster1, cluster2, float(toks[2]), float(toks[3]))
            else:
                add_clusters_to_dists(between_cluster_dist, cluster2, cluster1, float(toks[2]), float(toks[3]))

    cluster_ids = between_cluster_dist.keys()
    with open(os.path.join(out, "between_cluster_dist.csv"),"w") as f_out:
        f_out.write("cluster1, cluster2, core, core_max, accessory, accessory_max\n")
        for c1 in cluster_ids:
            for c2 in cluster_ids:
                if c2 <= c1:
                    continue
                core_max = max(between_cluster_dist[c1][c2]["core"])
                core_mean = mean(between_cluster_dist[c1][c2]["core"])
                acc_max = max(between_cluster_dist[c1][c2]["accessory"])
                acc_mean = mean(between_cluster_dist[c1][c2]["accessory"])
                f_out.write(",".join(map(str, [c1, c2, core_mean, core_max, acc_mean, acc_max])) + "\n")
    return


def run(args):
    args.out = os.path.abspath(args.out)
    args.metadata_file = os.path.abspath(args.metadata_file)
    fa_to_gff = assemblies_to_gffs(args.metadata_file)
    clusters, cluster_sizes = get_cluster_sizes(args.clusters_file, args.out)
    metadata_per_cluster(args.metadata_file, cluster_sizes, clusters, args.min_cluster_size, args.out)
    quit()
    within_cluster_dist = get_dist_within_cluster(clusters, cluster_sizes,
    args.dists_file, args.out, args.min_cluster_size, fa_to_gff)
    between_cluster_dist(clusters, within_cluster_dist, args.dists_file, fa_to_gff, args.out)
    return


def get_options():
    parser = argparse.ArgumentParser(description='Calculate the between and within cluster distances. Create files with roary jobs to run.')
    # input options
    parser.add_argument('--clusters_file', required=True,
                        type=str,
                        help='poppunk clusters output file')
    parser.add_argument('--dists_file',
                        required=True,
                        type=str,
                        help='output of extract_distances.py of core and accessory distances between every two strains')
    parser.add_argument('--metadata_file',
                        required=True,
                        type=str,
                        help='Metadata file containing the GFF files of the assemblies')
    parser.add_argument('--out', help='Name of output directory [%(default)s]', default = "dists_analysis",
                        type=str,
                        required=False)
    parser.add_argument('--min_cluster_size', help='Minimum size of a cluster to be used in further analysis [%(default)s]',
                        type=int,
                        default=30)
   	# parser.add_argument('--merge_files', action='store_true', help='Merge the distances and clusters file into one', default=False)
    return parser.parse_args()


if __name__ == "__main__":
    # get arguments from user
    start = time.time()
    options = get_options()
    run(options)
    end = time.time()
    print(end - start)
