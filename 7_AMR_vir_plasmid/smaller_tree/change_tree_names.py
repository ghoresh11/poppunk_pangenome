
d = {}
with open("smaller_tree.csv") as f:
	for line in f:
		toks = line.strip().split(",")
		filename = toks[3].split("/")[-1].replace(".gff","")
		d[filename] = toks[5]

with open("RAxML_bipartitions.tree") as f:
	for line in f:
		for filename in d:
			line = line.replace(filename, d[filename])
		print(line)