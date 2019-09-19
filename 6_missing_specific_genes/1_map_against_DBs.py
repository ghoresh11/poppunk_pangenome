import subprocess
import csv

## also need to run interpro scan

db_file = "DBs/ecosyc_genes.fa"
pathways_file = "DBs/MINIMAL_All_pathways_of_E._coli_K-12_substr._MG1655.txt"
results = {}

for gene_type in ["missing","specific"]:
    query_file = gene_type + ".fa"
    res_file = gene_type + "_blast_res.tab"

    # ## step 1: create a blast database
    subprocess.call(["makeblastdb", "-dbtype", "prot", "-in", db_file])

    ## step 2: blast against the dbfile
    subprocess.call(map(str, ["blastp", "-db", db_file,
                                   "-query", query_file, "-out", res_file,
                                   "-outfmt", "6 qseqid sseqid pident length qlen slen evalue bitscore", "-evalue", 0.01, "-num_threads",4]))

    ## step 3: see match gene to a process, add it to the existing "missing_genes_detailed_file"
    ## first create a dictionary of gene to process, then add to output
    ecosyc_dict = {}
    with open(pathways_file) as f:
        for toks in csv.reader(f, delimiter = "\t"):
            pathway = toks[0]
            genes = toks[1].replace("\"","").split(" // ")
            for g in genes:
                ecosyc_dict[g] = pathway

    with open(res_file) as f:  # add the edges to the graph
        for line in f:
            toks = line.strip().split()
            identity = float(toks[2])
            coverage = min(float(toks[3]) / float(toks[4]), float(toks[3]) / float(toks[5]))
            pangenome_gene = toks[0]
            pathway_gene = toks[1]
            if identity > 0.75 and coverage > 0.8:
                results[pangenome_gene] = ecosyc_dict[pathway_gene]

with open("ecosyc_results.csv","w") as out:
    out.write("Gene,Pathway\n")
    for g in results:
        out.write(g + "," + results[g] + "\n")
