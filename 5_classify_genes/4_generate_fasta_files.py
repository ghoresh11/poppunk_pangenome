from Bio.SeqIO.FastaIO import SimpleFastaParser

''' create a fasta file with all the sequences of a particular gene class'''

gene_to_sequence = {}
with open("/Users/gh11/poppunk_pangenome/4_pairwise_roary/231019_corrected/pan_genome_reference_prot.fa") as handle:
    for values in SimpleFastaParser(handle):
        gene_to_sequence[values[0].split()[0]] = values[1]


fasta_files = {}
with open("gene_classification.csv") as f:
    for line in f:
        toks = line.strip().split(",")
        gene_class = toks[1].replace(" ","_")
        gene_class = gene_class.replace("-","_")
        if gene_class not in fasta_files:
            fasta_files[gene_class] = open("sequences/" + gene_class + ".fa","w")
        fasta_files[gene_class].write(">" + toks[0] + "\n" + gene_to_sequence[toks[0]] + "\n")

for f in fasta_files:
    fasta_files[f].close()
