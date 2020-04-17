import os
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.Seq import translate

''' save a file "extract" with the member names of a cluster
to get all the sequences of a particular cluster'''

with open("leuO_extract") as f:
    for line in f:
        extract = line.strip().split()


cds_dir = "/nfs/pathogen004/gh11/40/"
files = os.listdir(cds_dir)
for f in files:
    if not f.endswith(".fa"):
        continue
    with open(os.path.join(cds_dir, f)) as handle:
        for values in SimpleFastaParser(handle):
            if values[0] in extract:
                print(">" + values[0] + "\n" + translate(values[1]))
