
print("## PANAROO")
cnt = 0
with open('no_dbscan/members.csv') as f:
	for line in f:
		if line.startswith("Gene"):
			continue
		toks = line.strip().split("\t")
		print(toks[0].split(",")[0])
		print(len(toks))

		cnt += 1
		if cnt == 20:
			break

print( "#### DBSCAN")
cnt = 0
with open('with_dbscan/members.csv') as f:
	for line in f:
		if line.startswith("Gene"):
			continue
		toks = line.strip().split("\t")
		print(toks[0].split(",")[0])
		print(len(toks))
		
		cnt += 1
		if cnt == 20:
			break
quit()
cnt = 0
print("## ROARY")
with open('/Users/gh11/poppunk_pangenome/4_pairwise_roary/180619/members.csv') as f:
	for line in f:
		if line.startswith("Gene"):
			continue
		toks = line.strip().split("\t")
		print(len(toks))
		cnt += 1
		if cnt == 10:
			break