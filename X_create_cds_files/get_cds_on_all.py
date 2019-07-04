import os
import subprocess
import sys

cnt = 0
with open("/lustre/scratch118/infgen/team216/gh11/e_coli_collections/FILTERED_MD_FINAL_ALL.tab") as f:
    for line in f:
        toks = line.strip().split("\t")
        if line.startswith("ID"):
            gff_loc = toks.index("New_annot_loc")
            poppunk_loc = toks.index("Poppunk_cluster")
            assembly_loc = toks.index("Assembly_Location")
            continue
        cluster = toks[poppunk_loc]
        if not os.path.exists(cluster):
            os.makedirs(cluster)

        gff_file = toks[gff_loc]
        curr_genome = gff_file.split("/")[-1].replace(".gff", "")
        assemblies = toks[assembly_loc].split(",")
        chosen_assembly = None
        for a in assemblies:
            curr = a.split("/")[-1].split(".")[0]
            if curr.lower() in gff_file.lower():
                chosen_assembly = a
        if chosen_assembly is None:
            sys.exit("gff: %s\nassemblies: %s" %(gff_file, str(assemblies)))


        mem = "100"
        lsf_prefix = ["bsub", "-q", "small", "-J", curr_genome, "-G", "team216","-o", curr_genome + ".o",
             "-e", curr_genome + ".e", '-R"select[mem>' + mem + '] rusage[mem='+ mem + ']"', '-M' + mem]

        command = map(str, lsf_prefix + ["python", "get_cds_files.py",
        "--gff_file", gff_file,
        "--fasta_file", chosen_assembly, "--out_dir", cluster])
        subprocess.call(command)
