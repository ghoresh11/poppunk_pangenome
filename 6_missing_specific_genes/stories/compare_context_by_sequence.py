import os
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.Seq import reverse_complement
import random

''' provide it with a list of genes, it will find the gene
upstream and downstream to the genes, this way I can see if the context is the same'''

def get_members_for_gene(gene):
	with open("/lustre/scratch118/infgen/team216/gh11/e_coli_collections/poppunk/new_roary/pairwise/analysis/071019_mode_rep/members.csv") as f:
		for line in f:
			toks = line.strip().split(",")
			if toks[0] == gene:
				members = toks[1].split("\t")
				break
	gene_per_cluster = {}
	for m in members:
		curr_cluster = m.split("_")[0]
		gene_name = "_".join(m.split("_")[1:])
		if curr_cluster not in gene_per_cluster:
			gene_per_cluster[curr_cluster] = []
		gene_per_cluster[curr_cluster].append(gene_name)

	members = {}
	roary_dir = "/lustre/scratch118/infgen/team216/gh11/e_coli_collections/poppunk/new_roary/"
	for curr_cluster in gene_per_cluster:
		members[curr_cluster] = []
		with open(os.path.join(roary_dir, curr_cluster, 'clustered_proteins')) as f:
			for line in f:
				toks = line.strip().split()
				curr_gene = toks[0].replace(":","")
				if curr_gene in gene_per_cluster[curr_cluster]:
					members[curr_cluster] += toks[1:]
	return members


def get_cluster_genomes(cluster):
	genomes = os.listdir(os.path.join("/nfs/pathogen004/gh11/",cluster))
	final_genomes = []
	for g in genomes:
		final_genomes.append(g.replace("_cds.fa",".gff"))
	random.shuffle(final_genomes)
	return final_genomes

def get_adjacent_regions(c, gene_group, out_up, out_down, genomes, members):
	cnt = 0
	for f in genomes:
		tmp = open("tmp.fa","w")
		contigs_to_identifiers = {}
		fasta = False
		with open(os.path.join("/nfs/pathogen004/gh11/gffs/",f)) as f_open:
			for line in f_open:
				if line.startswith("#"):
					continue
				if line.startswith(">"):
					fasta = True
				if fasta:
					tmp.write(line)

				toks = line.strip().split("\t")
				identifier = toks[-1].split(";")[0].replace("ID=","")
				name = toks[-1].split("Name=")[-1].split(";")[0]
				if len(name) == 0:
					name = toks[-1].split(";")[-1].replace("product=","")

				if identifier in members:
					if toks[0] not in contigs_to_identifiers:
						contigs_to_identifiers[toks[0]] = {}
					contigs_to_identifiers[toks[0]][identifier] = {"start": int(toks[3]), "stop": int(toks[4]), "strand": toks[5]}
		tmp.close()
		
		with open("tmp.fa") as handle:
			for values in SimpleFastaParser(handle):
				if values[0] not in contigs_to_identifiers:
					continue
				for gene in contigs_to_identifiers[values[0]]:
					seq = values[1]
					start = contigs_to_identifiers[values[0]][gene]["start"]
					stop = contigs_to_identifiers[values[0]][gene]["stop"]
					strand = contigs_to_identifiers[values[0]][gene]["strand"]
					pre = seq[max(0, start - 1000):start]
					post = seq[stop:min(stop + 1000, len(seq) )]
					identifier_out = ">" + c + "_" + gene_group + "(" + gene + ")" + "\n"

					cnt += 1
					if strand == "+":
						out_up.write(identifier_out + pre + "\n")
						out_down.write(identifier_out + post + "\n")
					else:
						out_up.write(identifier_out + reverse_complement(post) + "\n")
						out_down.write(identifier_out + reverse_complement(pre) + "\n")				
					
		if cnt >= 10:
			break



genes = ["flu_1****","group_1305","group_1559***","group_2980*******","group_349","group_3518************","group_40**","group_413**","group_5364****","group_5401***********","group_668**","group_6703***","group_77****","group_78***","group_8250**","group_84*****","group_869*","icsA_1**","icsA_2*","icsA_2****","icsA**","prn*","prn**","tibA**","tibA***","tibA*****","tibA******","yfaL_1**","yfaL**"]


name = "autotransporter"
out_up = open(name + "_upstream.fa","w")
out_down = open(name + "_downstream.fa","w")

for gene in genes:
	print("Getting context of gene %s..." %gene)
	members_per_cluster = get_members_for_gene(gene)

	for cluster in members_per_cluster:
		print("Getting context in cluster: %s..." %cluster)
		genomes = get_cluster_genomes(cluster)
		get_adjacent_regions(cluster, gene, out_up, out_down, genomes, members_per_cluster[cluster])
out_up.close()
out_down.close()