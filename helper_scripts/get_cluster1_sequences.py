import random
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.Seq import reverse_complement
import os

## I didn't manage to get the cluster 1 sequences for whatever reason
## SO i need to get a rep for each cluster from the list of genes that have it

### It's quite messy but should always work if this problem arrises again. The roary
## output is a bit unpredictable when parsing this is why this is the result

presence_absence_csv = "/lustre/scratch118/infgen/team216/gh11/e_coli_collections/poppunk/dists_analysis/roary_outputs/1_1550764508/gene_presence_absence.csv"
pan_genome_reference_out = "/lustre/scratch118/infgen/team216/gh11/e_coli_collections/poppunk/dists_analysis/roary_outputs/1_1550764508/pan_genome_reference.fa"

print("Choosing rep for each gene cluster...")
## Step 1: create a dictionary an isolate -> rep + gene name 
genes = {}
## genes -> genome -> name_in_gff -> name_in_presence_absence_file

with open(presence_absence_csv) as f:
	for line in f:
		if line.startswith("\"Gene"):
			continue
		toks = line.strip().split(",")
		gene_name = toks[0].replace("\"","")
		flag = False
		for i in range(len(toks)-1, 13, -1):
			r = toks[i]
			r = r.replace("\"","")
			r = r.split("\t")[-1]
			genome_name = r.split("_")[:-1]
			genome_name = "_".join(genome_name)
			if genome_name == "":
				continue
			if genome_name not in genes:
				genes[genome_name] = {}
			genes[genome_name][r] = gene_name
			flag = True
			break
		if not flag:
			print("Couldn't find rep for gene: %s" %toks[0])


def find_path(line, genes, genomes_to_loc):
	for genome in genes.keys():
		if genome in line:
			genomes_to_loc[genome] = line
			return 

print("Getting the GFF file locations...")
## Step 2: find the location of the GFF file of this isolate
genomes_to_loc = {}
with open("/lustre/scratch118/infgen/team216/gh11/e_coli_collections/poppunk/dists_analysis/gff_jobs/jobs_1.txt") as f:
	for line in f:
		line =  line.strip()
		find_path(line, genes, genomes_to_loc)
		


print("Getting the genes from the GFF files...")

## Step 3: for each isolate, find all the genes in the GFF file and save to a FASTA file
out = open(pan_genome_reference_out, "w")

for genome in genes:
	output_details = {}
	curr_genes = genes[genome].keys()
	with open(genomes_to_loc[genome]) as f:
		fasta = False
		tmp_out = open("fasta_tmp_out.fa","w")
		for line in f:
			if fasta:
				tmp_out.write(line)
				continue
			if line.startswith("##FASTA"):
				fasta = True
				continue
			if line.startswith("##"):
				continue
			toks = line.strip().split("\t")
			gene_name = toks[-1].split(";")[0]
			gene_name = gene_name.replace("ID=","")
			
			if gene_name not in curr_genes:
				continue
			if toks[0] not in output_details:
				output_details[toks[0]] = {}
			output_details[toks[0]][genes[genome][gene_name]] = {"strand":toks[6], "start": int(toks[3]), "stop": int(toks[4])}
		tmp_out.close()
		with open("fasta_tmp_out.fa") as handle:
			for values in SimpleFastaParser(handle):
				if values[0] in output_details:
					for gene in output_details[values[0]]:
						out.write(">" + gene + "\n")
						start = output_details[values[0]][gene]["start"]
						stop = output_details[values[0]][gene]["stop"]
						if output_details[values[0]][gene]["strand"] == "+":
							out.write(values[1][start-1:stop] + "\n")
						else:
							out.write(reverse_complement(values[1][start-1:stop]) + "\n")
out.close()
os.remove("fasta_tmp_out.fa")
