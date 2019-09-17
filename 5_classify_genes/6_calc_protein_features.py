import os
from Bio.Seq import translate
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.SeqUtils.ProtParam import ProteinAnalysis

''' instread of running this only on the rare variants, I will run it on all gene types and compare '''

def run(filename, fileout):
    cnt_stop_codon = 0
    cnt_ambiguous = 0
    out = open(fileout, "w")
    out.write("name,length,weight,aromaticity,instability,hydrophobicity\n")
    with open(filename) as handle:
        for values in SimpleFastaParser(handle):
            seq = values[1]
            if seq[-1] == "*":
                seq = seq[:-1]
            if "*" in seq:
                cnt_stop_codon += 1
                continue
            if "X" in seq:
                cnt_ambiguous += 1
                out.write(",".join(map(str,[values[0].split()[0],
                len(values[1])/3, "NA", "NA",
                "NA", "NA"])) + "\n")
                continue
            seq = ProteinAnalysis(seq)
            out.write(",".join(map(str,[values[0].split()[0],
            len(values[1])/3, seq.molecular_weight(), seq.aromaticity(),
            seq.instability_index(),
            seq.gravy()])) + "\n")
    out.close()
    print("Num premature stops: %d\nNum ambiguous: %d" %(cnt_stop_codon, cnt_ambiguous))
    return

if __name__ == "__main__":
    classification_dir = "/lustre/scratch118/infgen/team216/gh11/e_coli_collections/poppunk/classify_genes/prot/"
    all_gene_type_files = os.listdir(classification_dir)
    for gene_type_file in all_gene_type_files:
        if not gene_type_file.endswith(".fa"):
            continue
        print(gene_type_file + "...")
        curr_file_in = os.path.join(classification_dir, gene_type_file)
        curr_file_out = os.path.join("props",gene_type_file.replace("_prot.fa", "_props.csv"))
        run(curr_file_in, curr_file_out)
