import argparse
import os


def read_summary_files(input_dir):
    ''' read in the summary file so we have a dictionary
    with all the plasmidFinder identifiers and the members of that identfier'''
    summary_files = os.listdir(input_dir)
    plasmids = set()
    genomes = {}

    for f in summary_files:
        if not f.endswith("_summary.csv"):
            continue
        curr_cluster = f.split("_")[0]
        with open(os.path.join(input_dir, f)) as f_open:
            for line in f_open:
                toks = line.strip().split(",")
                if line.startswith("name"):
                    curr_plasmids = toks[1:]
                    for p in curr_plasmids:
                        plasmids.add(p)
                    continue
                name = toks[0].split("/")[-2]
                genomes[name] = []
                for i in range(1, len(toks)):
                    if toks[i] == "yes":
                        genomes[name].append(curr_plasmids[i-1])
    print(plasmids)
    return (plasmids, genomes)

def write_summary(plasmids, genomes):
    ''' write all to a file'''
    plasmids = list(plasmids)
    out = open("plasmid_summary_all.csv", "w")
    out.write("strain, cluster," + ",".join(plasmids) + "\n")
    for g in genomes:
        name = "_".join(g.split("_")[1:])
        cluster = g.split("_")[0]
        out.write(name + "," + cluster)
        for p in plasmids:
            if p in genomes[g]:
                out.write(",1")
            else:
                out.write(",0")
        out.write("\n")
    out.close()
    return


def run(args):
    plasmids, genomes = read_summary_files(os.path.abspath(args.input_dir))
    write_summary(plasmids, genomes)
    return


def get_options():
    parser = argparse.ArgumentParser(
        description='Summarise the ARIBA outputs')
    # input options
    parser.add_argument('--input_dir', required=False,
                        type=str, default="/lustre/scratch118/infgen/team216/gh11/e_coli_collections/plasmids/outputs/summary_files/",
                        help='path to ARIBA summary files [%(default)s]')

    return parser.parse_args()


if __name__ == "__main__":
    # get arguments from user
    options = get_options()
    run(options)
