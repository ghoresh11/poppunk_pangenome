import os
import subprocess

files = os.listdir(".")

for fasta_file in files:
    if not fasta_file.endswith(".fa"):
        continue    
    job_name = fasta_file.split(".")[0]
    queue = "parallel"
    mem = "8000"
    threads = "16"

    complete = []
    ## skip what had already completed
    if job_name in complete:
        continue

    lsf_prefix = ["bsub", "-q", queue, "-J", job_name, "-G", "team216", "-o", job_name + ".o",
                  "-e", job_name + ".e", '-R"select[mem>' + mem + '] rusage[mem=' + mem + '] span[hosts=1]"', '-M' + mem, "-n" + threads]

    command = ["emapper.py", "-d", "bact", "--dmnd_db", "/nfs/pathogen/eggnog/eggnog_proteins_v9.dmnd",
                "--data_dir","/nfs/pathogen/eggnog/",
               "--output", os.path.join("eggnog", job_name + ".txt"), "-m", "diamond", "-i", fasta_file, "--cpu", threads]
    subprocess.call(lsf_prefix + command)