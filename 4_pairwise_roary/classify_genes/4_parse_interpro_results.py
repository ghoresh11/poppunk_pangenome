import os
import re


## 1. constuct dictionary of
## interpro_scan -> gene name
## GO term -> gene names
## in the end make a table of gene_name -> GO terms (biological process, namespaces...) + interpro_scan (family, maybe domain, Homologous_superfamily)

class Entry:

    def __init__(self, name):
        self.name = name
        self.classification = "unknown"
        self.pfam_signature = set()
        self.biological_process = set()
        self.molecular_function = set()
        self.cellular_component = set()
        self.family = set()
        self.domain = set()
        self.superfamily = set()
        return

    def __str__(self):
        return "\t".join([self.name, self.classification, "/".join(list(self.pfam_signature)),
        "/".join(list(self.biological_process)), "/".join(list(self.molecular_function)),
        "/".join(list(self.cellular_component)), "/".join(list(self.superfamily)),
        "/".join(list(self.family)), "/".join(list(self.domain))]
        )

gff_files = os.listdir(".")
for gff_file in gff_files:
    if not gff_file.endswith("_result"):
        continue
    print(gff_file  + "...")
    fasta_file = os.path.join("/lustre/scratch118/infgen/team216/gh11/e_coli_collections/poppunk/classify_genes/prot", gff_file.replace("_result",""))
    print(fasta_file + "...")

    interpro_scan_ids = {   }
    go_terms = {}
    genes = {}

    name_to_gene = {}
    with open(fasta_file) as f:
        for line in f:
            if line.startswith(">"):
                toks = line.strip().split()
                entry_name = toks[0][1:]
                entry_name = re.sub(r"[^a-zA-Z0-9_\*-]+", '', entry_name) ## remove special characters
                genes[entry_name] = Entry(entry_name)

    with open(gff_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            if line.startswith(">"):
                break
            toks = line.strip().split("\t")
            entry_name = toks[0]
            if "InterPro" in line:
                ipr = line.split("IPR")[-1]
                ipr = ipr.split("\"")[0]
                ipr = "IPR" + ipr
                if ipr not in interpro_scan_ids:
                    interpro_scan_ids[ipr] = set()
                interpro_scan_ids[ipr].add(entry_name)
                ## get the interpro id and save in interpro scan ids
            if "signature_desc" in line: ## pfam
                ## get pfam desciptor and add to entry
                pfam = line.split("signature_desc=")[-1]
                pfam = pfam.split(";")[0]
                genes[entry_name].pfam_signature.add(pfam)
            if "Ontology_term" in line:
                gos = line.strip().split("Ontology_term=")[-1]
                gos = gos.split(";")[0]
                gos = gos.replace("\"","")
                gos = gos.split(",")
                for go in gos:
                    if go not in go_terms:
                        go_terms[go] = set()
                    go_terms[go].add(entry_name)



    ## step 2: go over the interpro scan file to extract the family etc.
    with open("050819_entry.list") as f:
        for line in f:
            toks = line.strip().split("\t")
            if toks[1] not in ["Homologous_superfamily","Family","Domain"]:
                continue
            if toks[0] not in interpro_scan_ids:
                continue
            for gene in interpro_scan_ids[toks[0]]:
                if toks[1] == "Homologous_superfamily":
                    genes[gene].superfamily.add(toks[2])
                elif toks[1] == "Family":
                    genes[gene].family.add(toks[2])
                else:
                    genes[gene].domain.add(toks[2])


    ## step 3: go over the GO term file and extact the namespaces
    flag = False
    with open("/lustre/scratch118/infgen/pathogen/pathpipe/go_ontology/gene_ontology_20180201.obo") as f:
        for line in f:
            if line.startswith("id:"):
                curr_term = line.strip().split()[-1]
                if curr_term not in go_terms:
                    flag = False
                else:
                    flag = True
                continue
            if flag and line.startswith("name:"):
                curr_name = line.strip().replace("name: ","")
                continue
            if flag and line.startswith("namespace"):
                curr_namespace = line.strip().replace("namespace: ","")
                if curr_namespace == "biological_process":
                    for gene in go_terms[curr_term]:
                        genes[gene].biological_process.add(curr_name)
                elif curr_namespace == "molecular_function":
                    for gene in go_terms[curr_term]:
                        genes[gene].molecular_function.add(curr_name)
                else:
                    for gene in go_terms[curr_term]:
                        genes[gene].cellular_component.add(curr_name)
                flag = False
    ## this are processess under the "biological process category which are useful for me"
    defined_processes = ["transmembrane transport",
    "DNA recombination", "DNA replication","methylation",
    "pathogenesis","cell adhesion","oxidation-reduction process", "regulation of transcription", "metabolic process",
    "pilus assembly", "DNA integration", "transposition", "proteolysis","pilus organization","cell motility","flagellum",
    "biosynthetic process","signal transduction","secretion","DNA modification","transport","repair","conjugation",
    "DNA restriction-modification system","DNA topological change","competence","CRISPR","biosynthetic process",
    "catabolic process","cell wall organization","phosphorylation","protein modification","DNA packaging","tRNA","viral"]

    defined_components = ["membrane","flagellum","pilus"]

    defined_words_in_family = ["secretion system","ROK family","Transposase", "Phage","enzyme","hydrolase",
    "secretion","biosynthesis","CRISPR","Fimbrial","capsid","permease","Biofilm","toxin","DNA-binding","relaxase",
    "secreted","Flagellar","tail","Plasmid","pilus","conjugative"]


    for gene in genes:
        gene = genes[gene]
        gene_line = str(gene).lower()
        if ("secretion system" in gene_line or "secretion-system" in gene_line):
            gene_line = gene_line.replace("-"," ")
            system_type = gene_line.split(" secretion system")[0]
            system_type = system_type.split("type ")[-1]
            gene.classification = "type " + system_type + " secretion system"
            continue
        if "phage" in gene_line or "viral" in gene_line or "virus" in gene_line or "capsid" in gene_line:
            gene.classification = "prophage"
            continue
        for p in defined_processes:
            if p.lower() in str(gene.biological_process).lower():
                 gene.classification = p
                 break ## breaks defined processes loop
        if gene.classification != "unknown":
            continue ## will break the bigger loop
        if len(gene.biological_process) > 0:
            gene.classification = list(gene.biological_process)[0]
            continue
        if len(gene.molecular_function) > 0:
            gene.classification = list(gene.molecular_function)[0]
            continue

        for c in defined_components:
            if c in str(gene.cellular_component).lower():
                gene.classification = c
                break
        if gene.classification != "unknown":
            continue
        for w in defined_words_in_family:
            if w.lower() in str(gene.family):
                gene.classification = w
                break
        if gene.classification != "unknown":
            continue

        if "plasmid" in gene_line:
            gene.classification = "plasmid"
        elif "conjuga" in gene_line:
            gene.classification = "conjugation"
        elif "mobilisation" in gene_line or "mobilization" in gene_line:
            gene.classification = "mobilisation"
        elif "toxin-antitoxin" in gene_line or "toxin antitoxin" in gene_line:
            gene.classification = "toxin-antitoxin system"
        elif "abc transporter" in gene_line:
            gene.classification = "transmembrane transport"
        elif "metabolic process" in gene_line:
            gene.classification = "metabolic process"
        elif gene.biological_process ==  "cell adhesion" and "fimbria" in gene_line:
            gene.biological = "cell adhesion"
        elif "transposase" in str(gene.pfam_signature).lower():
            gene.classification = "transposase"
        elif "resolvase" in str(gene.pfam_signature).lower():
            gene.classification = "resolvase"
        elif "membrane" in gene_line:
            gene.classification = "membrane"
        elif "crispr" in gene_line:
            gene.classification = "CRISPR"
        elif "toxin" in gene_line:
            gene.classification = "toxin"
        elif "pilus" in gene_line:
            gene.classification = "pilus"
        elif "duf" in gene_line:
            gene.classification = "DUF"
        elif sum([len(gene.pfam_signature), len(gene.biological_process), len(gene.molecular_function),
        len(gene.cellular_component), len(gene.family), len(gene.domain),len(gene.superfamily)]) > 0:
            gene.classification = "general function prediction"

    with open(gff_file.replace(".fa_result",".csv"),"w") as out:
        out.write("name\tclassification\tpfam_signature\tbiological_process\tmolecular_function\tcellular_component\tsuperfamily\tfamily\tdomain\n")
        for gene in genes:
            out.write(str(genes[gene]) + "\n")
