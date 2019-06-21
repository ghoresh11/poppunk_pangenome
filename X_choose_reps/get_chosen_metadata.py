import os

def get_input_dirs(input_dir):
    ''' check all directories in the input dir
    return a list of directories that have all the required files'''
    print("Getting the input directories...")
    input_dir = os.path.abspath(input_dir)
    directories = [d for d in os.listdir(input_dir)]
    dirs_to_return = {}
    for d in directories:
        d = os.path.join(input_dir, d)
        cluster = d.split("/")[-1].split("_")[0]
        if not os.path.isdir(d):
            continue
        # if cluster in ["51"]: ## my test rabbits
        #     continue
        if os.path.isfile(os.path.join(d,"tree_trimmed_list_X_10")):
            dirs_to_return[cluster] = os.path.join(d,"tree_trimmed_list_X_10")
    return dirs_to_return


dirs = get_input_dirs("/lustre/scratch118/infgen/team216/gh11/e_coli_collections/poppunk/new_roary/")

chosen = {}
for d in dirs:
	with open(dirs[d]) as f:
		for line in f:
			name = line.strip().lower().replace("_trimmed","").replace("_contigs_pacbio","").replace("_contigs_canu_1_6","").replace("_contigs_hgap_4_0","")
			chosen[name] = d.split("/")[-1].split("_")[0]

minimized_metadata = open("210619_chosen_treemer.csv","w")
minimized_metadata.write(",".join(["Name","Assembly", "Annotation", "Reads","popppunk_cluster","ST","Year","Pathotype","Country","Continent","Isolation","Publication", "num_contigs", "length"]) + "\n")


with open("/lustre/scratch118/infgen/team216/gh11/e_coli_collections/FILTERED_MD_FINAL_ALL.tab") as f:
	for line in f:
		toks = line.strip().split("\t")
		if line.startswith("ID"):
			assembly_loc = toks.index("Assembly_Location")
			gff_loc = toks.index("New_annot_loc")
			reads_loc = toks.index("Reads_Location")
			st_loc =  toks.index("ST")
			pathotype_loc = toks.index("Pathotype")
			country_loc = toks.index("Country")
			continent_loc = toks.index("Continent")
			isolation_loc = toks.index("Isolation")
			publication_loc = toks.index("Publication")
			num_contigs_loc = toks.index("num_contigs")
			length_loc = toks.index("length")
			year_loc = toks.index("Year")
			popppunk_cluster = toks.index("Poppunk_cluster")
			run_id = toks.index("Run_ID")
			continue
		name = None
		#print(toks[run_id].lower())
		

		if toks[0].lower() in chosen:
			name = toks[0].lower()
		elif toks[run_id].lower() in chosen:
			name = toks[run_id].lower()
		elif toks[1].lower() in chosen:
			name = toks[1].lower()
		else:
			continue

		minimized_metadata.write(",".join([toks[0], toks[assembly_loc].split(",")[0], toks[gff_loc], toks[reads_loc].replace(",",";"), chosen[name], toks[st_loc], toks[year_loc],
		toks[pathotype_loc], toks[country_loc], toks[continent_loc], toks[isolation_loc].replace(",","-"), toks[publication_loc].replace(",","-"), toks[num_contigs_loc], toks[length_loc]]) + "\n")
		chosen[name] = "DONE"

minimized_metadata.close()

for c in chosen:
	if chosen[c] != "DONE":
		print(c)