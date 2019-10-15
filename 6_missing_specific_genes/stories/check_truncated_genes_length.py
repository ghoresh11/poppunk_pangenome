import sys

''' requires the fasta file of a gene cluster taken from the NFS space
go over all the sequences to see what their length is'''

trunc_gene = sys.argv[1]

counts = {}
total = 0
with open(trunc_gene) as f:
	for line in f:
		if line.startswith(">"):
			continue
		curr_length = len(line.strip())/3
		if curr_length not in counts:
			counts[curr_length] = 0
		counts[curr_length] += 1
		total += 1

for c in counts:
	counts[c] = counts[c]/float(total)
print(counts)