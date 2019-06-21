
cd ${1}

#snp_sites core_gene_alignment.aln -m -o core_gene_alignment_snps.aln

#FastTree -nt core_gene_alignment_snps.aln > tree

## This will run treemer and keep only 10 leaves
python ~/Treemmer/Treemmer.py -np -X 10 tree