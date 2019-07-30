import os
import subprocess

# Input/output:
#   -i INPUT_FILES [INPUT_FILES ...], --input INPUT_FILES [INPUT_FILES ...]
#                         input GFF3 files (usually output from running Prokka)
#   -o OUTPUT_DIR, --out_dir OUTPUT_DIR
#                         location of an output directory

# Matching:
#   -c ID, --threshold ID
#                         sequence identity threshold (default=0.95)
#   -f FAMILY_THRESHOLD, --family_threshold FAMILY_THRESHOLD
#                         protein family sequence identity threshold
#                         (default=0.7)
#   --len_dif_percent LEN_DIF_PERCENT
#                         length difference cutoff (default=0.95)
#   --merge_paralogs      don't split paralogs

# Graph correction:
#   --mode {strict,moderate,relaxed}
#                         the stringency mode at which to run panaroo. One of
#                         'strict', 'moderate' or 'relaxed' (default='strict')
#   --min_trailing_support MIN_TRAILING_SUPPORT
#                         minimum cluster size to keep a gene called at the end
#                         of a contig
#   --trailing_recursive TRAILING_RECURSIVE
#                         number of times to perform recursive triming of low
#                         support nodes near the end of contigs
#   --edge_support_diff EDGE_SUPPORT_DIFF
#                         maximum fraction difference between an edge's support
#                         and those of the nodes it connects
#   --remove_by_consensus
#                         if a gene is called in the same region with similar
#                         sequence a minority of the time, remove it
#   --high_var_flag CYCLE_THRESHOLD_MIN
#                         minimum number of nested cycles to call a highly
#                         variable gene region (default = 5).
#   --min_edge_support_sv MIN_EDGE_SUPPORT_SV
#                         minimum edge support required to call structural
#                         variants in the presence/absence sv file

# Gene alignment:
#   -a ALN, --alignment ALN
#                         Output alignments of core genes or all genes. Options
#                         are 'core' and 'pan'. Default: 'None'
#   --aligner ALR         Specify an aligner. Options:'prank', 'clustal', and
#                         default: 'mafft'
#   --core_threshold CORE
#                         Core-genome sample threshold (default=0.95)


jobs_dir= "/lustre/scratch118/infgen/team216/gh11/e_coli_collections/poppunk/new_roary/new_jobs_corrected/"
panaroo_out_dir= "/lustre/scratch118/infgen/team216/gh11/e_coli_collections/poppunk/panaroo/"
jobs = os.listdir(jobs_dir)

complete = ["30", "8", "43", "40", "21"]
failed = ["1"]
cnt = 0
for jobs_file in jobs:
    if not jobs_file.endswith(".txt"):
        continue

    cluster_id = jobs_file.split("_")[1].replace(".txt","")

    if cluster_id in complete:
        continue

    if cluster_id not in failed:
        continue

    gffs = []
    with open(os.path.join(jobs_dir, jobs_file)) as f:
        for line in f:
            gffs.append(line.strip())

    queue = "normal"
    mem = 2500

    if int(cluster_id) < 7:
        queue = "long"
        mem = 3000

    if int(cluster_id) == 1:
        mem = 12000
        queue = "basement"

    mem = str(mem)
    ## create output directory for this cluster
    outdir = os.path.join(panaroo_out_dir, cluster_id)
    try:
        os.makedirs(outdir)
    except Exception:
        pass

    job_name = "cluster_"+cluster_id

    lsf_prefix = ["bsub", "-q", queue, "-J", job_name, "-G", "team216","-o", job_name + ".o",
     "-e", job_name + ".e", '-R"select[mem>' + mem + '] rusage[mem='+ mem + ']"', '-M' + mem]

    command = ["panaroo","-i"] + gffs + ["-o", cluster_id, "--mode", "relaxed"]
    ## submit the job
    subprocess.call(lsf_prefix + command)
