
cd ${1}

#snp_sites core_gene_alignment.aln -m -o core_gene_alignment_snps.aln

FastTreeMP -mlnni 4 -nt core_gene_alignment_snps.aln > tree

## This will run treemer and keep only 10 leaves
runjob5 -e more_cpus.e -o more_cpus.o -n6 -R"span[hosts=1]" python ~/Treemmer/Treemmer.py --cpu 6 -np -X 10 tree