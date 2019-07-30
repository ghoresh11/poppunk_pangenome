import os
import gffutils
import sys
import subprocess

def get_assembly_files(gff_files, md_file="/lustre/scratch118/infgen/team216/gh11/e_coli_collections/FILTERED_MD_FINAL_ALL.tab"):
    gff_to_fasta = {}
    with open(md_file) as f:
        for line in f:
            toks = line.strip().split("\t")
            if line.startswith("ID"):
                gff_loc = toks.index("New_annot_loc")
                poppunk_loc = toks.index("Poppunk_cluster")
                assembly_loc = toks.index("Assembly_Location")
                continue

            gff_file = toks[gff_loc]
            if gff_file not in gff_files:
                continue
            assemblies = toks[assembly_loc].split(",")
            chosen_assembly = None
            for a in assemblies:
                curr = a.split("/")[-1].split(".")[0]
                if curr.lower() in gff_file.lower():
                    chosen_assembly = a
            gff_to_fasta[toks[gff_loc]] = chosen_assembly
    return gff_to_fasta

def parse_cdhit_clusters(in_file):
    member_to_rep = {}
    curr_members = []
    curr_rep = None
    with open(in_file) as f:
        for line in f:
            if line.startswith(">"):
                for m in curr_members:
                    member_to_rep[m] = curr_rep
                curr_members = []
                curr_rep = None
                continue
            member = line.strip().split("...")[0].split(">")[-1]
            curr_members.append(member)
            if "*" in line:
                curr_rep = member
    ## don't forget the last one
    for m in curr_members:
        member_to_rep[m] = curr_rep
    return member_to_rep

def run_blast(fasta_file):
    ''' run all against all blast on input file'''
    print("Running BLAST...")
    cdhit_tmp_out = "cdhit_tmp_out"
    subprocess.call(["cd-hit-est", "-i", fasta_file ,"-o", "cdhit_tmp_out", "-c", "1", "-n","10", "-T", "4", "-d", "0","-A", "1"])
    member_to_rep = parse_cdhit_clusters(cdhit_tmp_out + ".clstr")
    makeblastdb_command = map(
        str, ["makeblastdb", "-in", cdhit_tmp_out, "-dbtype", "nucl"])
    blastn_command = map(str, ["blastn", "-db", cdhit_tmp_out,
                               "-query", cdhit_tmp_out, "-out",
                              "tmp_blast_results.tab",
                               "-outfmt", "6 qseqid sseqid pident length qlen slen evalue bitscore", "-evalue","1", "-num_threads", "4"])
    subprocess.call(makeblastdb_command)
    subprocess.call(blastn_command)
    return member_to_rep

def get_gene_lengths(ids_of_interest, filename = "/lustre/scratch118/infgen/team216/gh11/e_coli_collections/poppunk/new_roary/new_jobs_corrected/jobs_30.txt"):
    lengths = {}
    seqs = {}
    gff_files =  [line.strip('\n') for line in open(filename)]
    gff_to_fasta = get_assembly_files(gff_files)
    tmp_out = open("tmp_out.fa","w")
    for gff_file in gff_files:
        if not gff_file.endswith(".gff"):
            continue
        print("reading file: %s" %gff_file)
        db = gffutils.create_db(gff_file, dbfn=  gff_file + '.db', force=True)
        for i in db.features_of_type('CDS'):
            if db[i].id not in ids_of_interest:
                continue
            lengths[db[i].id] = db[i].stop - db[i].start
            seq = db[db[i].id].sequence(gff_to_fasta[gff_file])
            tmp_out.write(">" + db[i].id + "\n" + seq + "\n")
        os.remove(gff_file + ".db")
    tmp_out.close()
    member_to_rep = run_blast("tmp_out.fa")
    res = {}
    with open("tmp_blast_results.tab") as f:
        for line in f:
            toks = line.strip().split()
            if toks[0] not in res:
                res[toks[0]] = {}
            res[toks[0]][toks[1]] = toks[2]
    return lengths, res, member_to_rep



keys = ["1000", "0111"]
special_edges = {}
for k in keys:
    special_edges[k] = []

with open("edgelist") as f:
    for line in f:
        toks = line.strip().split()
        val = "".join(toks[2:])
        if val in keys:
            special_edges[val].append(toks[:2])

flatten = lambda l: [item for sublist in l for item in sublist]
ids_of_interest = set(flatten(flatten(special_edges.values())))
lengths, blast_res, member_to_rep = get_gene_lengths(ids_of_interest)


out = open("summary_diff_edges_updated.csv", "w")
out.write("edge_type,refound,length_diff, alignment_score\n")
for k in keys:
    for edge in special_edges[k]:
        if edge[0] not in lengths or edge[1] not in lengths:
            if "refound" in edge[0] or "refound" in edge[1]:
                out.write(k + ",1,NA,NA\n")
            else:
                sys.exit("ERROR: could not find edge length %s in GFF files" %str(edge))
        else:
            diff = abs(lengths[edge[0]] - lengths[edge[1]])
            rep1 = member_to_rep[edge[0]]
            rep2 = member_to_rep[edge[1]]
            if rep1 is None or rep2 is None:
                print("first: %s (%s), second: %s (%s)" %(edge[0],  str(rep1), edge[1], str(rep2)))
            elif rep1 == rep2:
                    ident = 100
            elif rep1 in blast_res and rep2 in blast_res[rep1]:
                ident = blast_res[rep1][rep2]
            elif rep2 in blast_res and rep1 in blast_res[rep2]:
                ident = blast_res[rep2][rep1]
            else:
                ident = "0"
                continue
            out.write(k + ",0," + str(diff) + "," + str(ident) + "\n")
out.close()
