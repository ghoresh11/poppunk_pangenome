import os

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
	return final_genomes


def get_adjacent_genes(c, out, genomes, members):
	cnt = 0
	for f in genomes:
		double_prev_id = ""
		prev_id = ""
		curr_id = ""
		double_prev_name, prev_name, curr_name = "","",""
		with open(os.path.join("/nfs/pathogen004/gh11/gffs/",f)) as f_open:
			for line in f_open:
				if line.startswith("#"):
					continue
				if line.startswith(">"):
					break

				toks = line.strip().split("\t")
				identifier = toks[-1].split(";")[0].replace("ID=","")
				name = toks[-1].split("Name=")[-1].split(";")[0]
				if len(name) == 0:
					name = toks[-1].split(";")[-1].replace("product=","")

				drouble_prev_id = prev_id
				prev_id = curr_id
				curr_id = identifier

				double_prev_name = prev_name
				prev_name = curr_name
				curr_name = name

				if prev_id in members:
					if toks[6] == "+":
						out.write(c + "," + f + "," + double_prev_name + "," + prev_id + "," + curr_name +  "\n")
					else:
						out.write(c + "," + f + "," + curr_name + "," + prev_id + "," + double_prev_name  + "\n")
					cnt += 1
		if cnt >= 10:
			break


geneA = "yehB_2****"
geneB = "yehB_1**"
geneC = "group_1898********"

geneA = "aaaT"
for gene in [geneA]:
	print("Getting context of gene %s..." %gene)
	members_per_cluster = get_members_for_gene(gene)
	out = open(gene + ".csv","w")
	out.write("Genome,Cluster,Before,Gene,After\n")
	for cluster in members_per_cluster:
		print("Getting context in cluster: %s..." %cluster)
		genomes = get_cluster_genomes(cluster)
		get_adjacent_genes(cluster, out, genomes,members_per_cluster[cluster])
	out.close()


